#!/usr/bin/env python3
"""isoforge annotate — transcript models, their ends, counts and fusions from long reads.

    isoforge annotate --reads a.fastq.gz b.fastq.gz --genome hg38.fa --annotation gencode.gtf --o results/project
    isoforge annotate --bam a.bam b.bam --genome hg38.fa --annotation gencode.gtf --o results/project

Read files are aligned first; sorted BAM files are read as they are. The run builds and admits the transcript
models and writes the summary, the filter list, the abundance table, the GTF, the read table and the fusion calls.
"""


import multiprocessing
import os
import sys
import time
from pathlib import Path

from isoforge import options
from isoforge.genome_fasta import ensure_fai
from isoforge.isoforge_align import align_samples, is_read_file, remove_alignments, sample_name
from isoforge.isoforge_init import index_for_gtf, index_metadata
from isoforge.options import ADMISSION_CLASSES, SIBLING_EDIT_RATIOS

__version__ = "1.0.0"






def _admission_kwargs(args):
    """The learned-admission settings for isoforge_core.post_filter."""
    model = getattr(args, "admission_model", None) or str(Path(__file__).resolve().parent / "data" / "admission_v1.json")
    if not Path(model).exists():
        raise SystemExit(f"admission model not found: {model}")
    guard = [float(x) for x in str(getattr(args, "admission_guard", "0.6,0.6,0.1,0.8,0.6")).split(",")]
    if not 3 <= len(guard) <= 5:
        raise SystemExit("--admission_guard takes TAILED,TES,IP[,TSS_OUTSIDE[,TSS]]")
    guard += [0.8, 0.6][len(guard) - 3:]
    tailed, tes, ip, tss, tss_all = guard
    by_class = []
    spec = str(getattr(args, "admission_classes", ADMISSION_CLASSES) or "none")
    if spec.lower() != "none":
        for item in spec.split(","):
            part = item.strip().split(":")
            try:
                # 'never': no candidate of the class is added (ADD), no model of it removed (VETO)
                add = float("inf") if part[1].lower() == "never" else float(part[1])
                veto = 0.0 if part[2].lower() == "never" else float(part[2])
                by_class.append((part[0], add, veto))
            except (IndexError, ValueError):
                raise SystemExit(f"--admission_classes takes CLASS:ADD:VETO[,...] or 'none'; got {item!r}")
    return dict(admission_model=model, admission_add=float(getattr(args, "admission_add", 0.6)),
                admission_drop=float(getattr(args, "admission_drop", 0.5)), admission_by_class=by_class, admission_guard_tailed=tailed,
                admission_guard_tes=tes, admission_guard_ip=ip, admission_guard_tss=tss, admission_guard_tss_all=tss_all,
                admission_min_datasets=int(getattr(args, "admission_min_datasets", 2)))

def _sibling_ratio_kwargs(args):
    """The read ratios of the four sibling edits for isoforge_core.post_filter (--sibling_edit_ratios)."""
    spec = str(getattr(args, "sibling_edit_ratios", SIBLING_EDIT_RATIOS) or SIBLING_EDIT_RATIOS)
    try:
        terminal, microexon, shift, intron = (float(x) for x in spec.split(","))
    except ValueError:
        raise SystemExit(f"--sibling_edit_ratios takes TERMINAL,MICROEXON,SHIFT,INTRON; got {spec!r}")
    return dict(sibling_ratio_terminal=terminal, sibling_ratio_microexon=microexon, sibling_ratio_shift=shift,
                sibling_ratio_intron=intron)


def _fusion_rescue_bams(args, database_dir=None):
    """(dataset, BAM) pairs for the fusion rescue, the run's sorted and indexed BAMs: from --fusion_rescue_config,
    else the BAMs the run read (args.run_bams), else the list `isoforge annotate` left in the run database
    (bam_paths.csv). None when the rescue is off or nothing can be read; `notes` then says why."""
    if not getattr(args, "fusion_rescue", True):
        return None, None
    if not getattr(args, "genome_fasta", None):
        return None, "no genome FASTA (--genome)"
    pairs = []
    sheet = getattr(args, "fusion_rescue_config", None)
    if sheet and os.path.exists(sheet):
        with open(sheet) as fh:
            for line in fh:
                f = [x.strip() for x in line.rstrip("\n").split(",")]
                if len(f) >= 4 and f[0] and not f[0].startswith("#"):
                    pairs.append((f[0], f[3]))
    elif getattr(args, "run_bams", None):
        pairs.extend(args.run_bams)
    elif database_dir is not None and (Path(database_dir) / "bam_paths.csv").exists():
        with open(Path(database_dir) / "bam_paths.csv") as fh:
            for line in fh:
                f = line.rstrip("\n").split(",", 1)
                if len(f) == 2 and f[0]:
                    pairs.append((f[0], f[1]))
    if not pairs:
        return None, "the run's BAMs are not known (pass --fusion_rescue_config)"
    indexed = [(d, b) for d, b in pairs if os.path.exists(b) and any(os.path.exists(b + x) for x in (".bai", ".csi"))]
    if not indexed:
        return None, "the run's BAMs are missing or not indexed"
    note = None if len(indexed) == len(pairs) else f"{len(pairs) - len(indexed)} of {len(pairs)} BAMs missing or not indexed"
    return indexed, note


ANNOTATE_EPILOG = """\
Examples:
  isoforge annotate --reads a.fastq.gz b.fastq.gz --genome hg38.fa --annotation gencode.v49.gtf \\
      --o results/project
  isoforge annotate --bam a.bam b.bam --genome hg38.fa --annotation gencode.v49.gtf --o results/project

Sample sheet (--config), one sample per line, with a read file or a BAM:
  Control_rep1,Control replicate 1,ONT,/data/control_rep1.fastq.gz

`isoforge annotate --help-all` lists every option.
"""


def get_args():
    parser = options.parser_for(
        "annotate",
        usage="isoforge annotate (--reads FASTQ [FASTQ ...] | --bam BAM [BAM ...] | --config CSV)\n"
              "                         --genome FASTA (--annotation GTF | --db DIR) --o PREFIX [options]",
        description="Transcript models with their start and poly(A) sites, read counts, small RNAs and gene fusions\n"
                    "from long reads, given as FASTQ files or as sorted BAM files.",
        epilog=ANNOTATE_EPILOG,
        show_all="--help-all" in sys.argv,
    )
    args = parser.parse_args()
    given = [name for name, value in (("--reads", args.reads), ("--bam", args.bam), ("--config", args.config_file))
             if value]
    if len(given) != 1:
        parser.error("give the reads one way: --reads, --bam or --config"
                     + (f" (got {' and '.join(given)})" if given else ""))
    files = args.reads or args.bam or []
    if args.names and len(args.names) != len(files):
        parser.error("--names needs one name per file of --reads or --bam")
    if args.dataset_name and (args.names or len(args.bam or []) != 1):
        parser.error("--dataset_name names a single --bam; use --names")
    options.take_genome_fasta(args)
    options.name_outputs(args, ("summarize", "filter", "abundance", "gtf", "read_annot", "all_post"))
    if args.sj_correction_window is not None and args.sj_correction_window < 0:
        parser.error("--sj_correction_window must be nonnegative")
    if args.min_count_supported != 0 and args.min_count_supported < 5:
        parser.error("--minCountSupported must be 0 (disabled) or at least 5")
    if not args.reference_end_support:
        # the two refinements of annotated-end evidence follow it
        args.known_singleton_ends = False
        args.ism_intact_end_support = False
    if args.chain_first and args.approx_chain_tolerance <= 0:
        parser.error("--chain_first requires --approx_chain_tolerance > 0")
    if args.assign_unique_ism == "off":
        args.ism_intact_end_support = False
    if not args.evidence_tss_rescue:
        args.ism_evidence_rescue = False
    if args.approx_chain_no_lock and args.approx_chain_tolerance <= 0:
        parser.error("--approx_chain_no_lock requires --approx_chain_tolerance > 0")
    if 0 < args.approx_chain_tolerance <= args.max_extended_snap:
        parser.error("--approx_chain_tolerance must exceed --max_extended_snap "
                     f"({args.max_extended_snap}); below it, vertex snapping has "
                     "already absorbed the offset and the chain test adds nothing")
    return args


def discover_sj_files(sj_dir: str) -> list:
    """Find all STAR SJ.out.tab files in a directory."""
    sj_path = Path(sj_dir)
    if not sj_path.is_dir():
        print(f"ERROR: --sj_dir is not a directory: {sj_dir}", file=sys.stderr)
        sys.exit(1)

    files = sorted(sj_path.glob("*SJ.out.tab"))
    if not files:
        print(f"WARNING: No *SJ.out.tab files found in {sj_dir}", file=sys.stderr)
    return [str(f) for f in files]


def parse_config_file(config_file: str) -> list:
    """Parse config CSV → list of (file, dataset_name, description, platform); the file is a BAM or a read file."""
    if not Path(config_file).exists():
        print(f"ERROR: Config file not found: {config_file}", file=sys.stderr)
        sys.exit(1)

    entries = []
    with open(config_file) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) == 4:
                name, desc, platform, bam = parts
            elif len(parts) == 2:
                name, bam = parts
                desc, platform = "", "ONT"
            else:
                print(
                    f"ERROR: Config line {lineno}: expected 2 or 4 fields, got {len(parts)}: {line}",
                    file=sys.stderr,
                )
                sys.exit(1)

            if not Path(bam).exists():
                print(f"ERROR: Config line {lineno}: file not found: {bam}", file=sys.stderr)
                sys.exit(1)
            entries.append((bam, name, desc, platform))

    if not entries:
        print("ERROR: No valid entries in config file!", file=sys.stderr)
        sys.exit(1)

    return entries


def input_samples(args) -> list:
    """(file, sample name) per sample, from --reads, --bam or the sample sheet."""
    if args.config_file:
        return [(path, name) for path, name, _, _ in parse_config_file(args.config_file)]
    files = args.reads or args.bam
    for path in files:
        if not Path(path).exists():
            print(f"ERROR: file not found: {path}", file=sys.stderr)
            sys.exit(1)
    if args.names:
        names = args.names
    elif args.dataset_name:
        names = [args.dataset_name]
    else:
        names = [sample_name(path) for path in files]
    return list(zip(files, names))


def run_post_processing(annot_dir, database_dir, out_prefix, args, project_name):
    """Run selected post-processing steps using the Rust engine.
    
    Args:
        annot_dir: Path to annotation Parquet directory
        database_dir: Path to run-time Parquet database directory
        out_prefix: Full output prefix path (e.g. /results/project)
        args: Parsed CLI args
        project_name: Project name (basename of prefix)
    """
    import isoforge_core

    # Load PostDB once in Rust — all data shared across steps
    print("\n[Post] Loading PostDB from Rust engine...")
    t0 = time.time()
    db_handle = isoforge_core.post_load_db(str(annot_dir), str(database_dir))
    hyper_fits = os.environ.get("ISOFORGE_HYPER_FITS")
    if hyper_fits:
        t0 = time.time()
        n_reads, n_fits, n_with = isoforge_core.post_hyper_fits(db_handle, hyper_fits)
        print(f"[Post] Read hypergraph: {n_fits:,} fits of {n_with:,} of {n_reads:,} read classes -> {hyper_fits} ({time.time()-t0:.1f}s)")
    print(f"       PostDB loaded in {time.time()-t0:.1f}s")
    # Annotated forms are split off unless annotated records are written at their annotated ends.
    if not getattr(args, "known_ends_as_annotated", False):
        n_split = isoforge_core.post_split_annotated_forms(
            db_handle,
            known_end_snap=getattr(args, "known_end_snap", 0),
            known_end_max_extension=getattr(args, "known_end_max_extension", 0),
            known_end_extension_min_fraction=getattr(args, "known_end_extension_min_fraction", 0.5),
        )
        print(f"       Annotated forms kept beside read-derived ends: {n_split} records split")

    try:
        if args.summarize or args.all_post:
            print("\n[Post] Summarizing datasets...")
            t0 = time.time()
            isoforge_core.post_summarize(db_handle, str(out_prefix))
            print(f"       Done in {time.time()-t0:.1f}s")

        # ── Filtering (always runs when abundance/gtf/all requested) ──────
        whitelist_file = None
        need_filter = args.filter or args.abundance or args.gtf or args.read_annot or args.all_post
        if need_filter:
            print("[Post] Filtering transcripts...")
            t0 = time.time()
            whitelist_file = f"{out_prefix}_filter.csv"
            isoforge_core.post_filter(
                db_handle, str(whitelist_file),
                max_frac_a=args.max_frac_a_filter,
                min_count=args.min_count,
                min_datasets=args.min_datasets,
                allow_genomic=args.allow_genomic,
                filter_known=args.filter_known,
                min_count_known=args.min_count_known,
                min_tss_support=args.min_tss_support,
                tss_upstream_rescue=getattr(args, "tss_upstream_rescue", False),
                evidence_tss_rescue=args.evidence_tss_rescue,
                junction_genome=args.genome_fasta if args.validate_novel_junctions else None,
                junction_scores=getattr(args, "splice_scores", None) if args.validate_novel_junctions else None,
                junction_sj_files=(args.sj_files or []) + (discover_sj_files(args.sj_dir) if args.sj_dir else []),
                min_count_supported=args.min_count_supported,
                require_phased_novel_ends=args.require_phased_novel_ends,
                phased_ends_spliced=getattr(args, "phased_ends_spliced", False),
                junction_prune_radius=args.junction_prune_radius,
                junction_prune_ratio=args.junction_prune_ratio,
                min_novel_locus_fraction=args.min_novel_locus_fraction,
                frac_a_annotated_end_window=args.frac_a_annotated_end_window,
                novel_exact_support=getattr(args, "novel_exact_support", 0),
                novel_exact_tol=getattr(args, "novel_exact_tol", 0),
                ism_evidence_rescue=getattr(args, "ism_evidence_rescue", False),
                ism_rescue_min_reads=getattr(args, "ism_rescue_min_reads", 5),
                ism_rescue_min_distance=getattr(args, "ism_rescue_min_distance", 0),
                min_count_nic=getattr(args, "min_count_nic", 0),
                weak_source_tier=getattr(args, "weak_source_tier", False),
                sub_path_min_fraction=getattr(args, "sub_path_min_fraction", 0.0),
                sub_path_tail_evidence=getattr(args, "sub_path_tail_evidence", False),
                sibling_edits=getattr(args, "sibling_edits", True),
                **_sibling_ratio_kwargs(args),
                known_end_evidence=getattr(args, "known_end_evidence", 0),
                known_end_evidence_window=getattr(args, "known_end_evidence_window", 50),
                min_count_end_variant=getattr(args, "min_count_end_variant", 0),
                end_variant_min_datasets=getattr(args, "end_variant_min_datasets", 2),
                cap_tss_min_fraction=getattr(args, "cap_tss_min_fraction", 0.15),
                cap_tss_min_reads=getattr(args, "cap_tss_min_reads", 3),
                quant_requires_end_support=getattr(args, "quant_requires_end_support", False),
                start_model=getattr(args, "start_model", False),
                start_min_llr=getattr(args, "start_min_llr", 20.0),
                known_phased_reads=getattr(args, "known_phased_reads", 0),
                pool_datasets=getattr(args, "pool_datasets", False),
                **_admission_kwargs(args),
            )
            n_back = isoforge_core.post_reclaim_annotated_records(db_handle, str(whitelist_file))
            print(f"       Annotated records reported in place of their only passing end variant: {n_back}")
            print(f"       Done in {time.time()-t0:.1f}s")

        if args.abundance or args.all_post:
            print("[Post] Generating abundance matrix...")
            t0 = time.time()
            if whitelist_file and Path(whitelist_file).exists():
                isoforge_core.post_abundance(
                    db_handle, str(out_prefix), args.genome,
                    whitelist_file=str(whitelist_file),
                    quant_requires_end_support=getattr(args, "quant_requires_end_support", False),
                    max_frac_a=args.max_frac_a_filter,
                    frac_a_annotated_end_window=args.frac_a_annotated_end_window,
                    quant_reassign=True,
                    quant_partial_reads=False,
                    quant_hypergraph=True,
                )
            if getattr(args, 'unfiltered', False):
                isoforge_core.post_abundance(db_handle, str(out_prefix), args.genome,
                                             quant_requires_end_support=getattr(args, "quant_requires_end_support", False),
                                             max_frac_a=args.max_frac_a_filter,
                                             frac_a_annotated_end_window=args.frac_a_annotated_end_window,
                                             quant_reassign=True,
                                             quant_partial_reads=False,
                                             quant_hypergraph=False)
            print(f"       Done in {time.time()-t0:.1f}s")

        if args.gtf or args.all_post:
            print("[Post] Generating GTF (filtered transcripts only)...")
            t0 = time.time()
            isoforge_core.post_gtf(
                db_handle, str(out_prefix), args.genome,
                whitelist_file=str(whitelist_file) if whitelist_file else None,
                observed_only=True,
                known_end_snap=getattr(args, "known_end_snap", 0),
                known_end_max_extension=getattr(args, "known_end_max_extension", 0),
                known_end_extension_min_fraction=getattr(args, "known_end_extension_min_fraction", 0.5),
                known_ends_as_annotated=getattr(args, "known_ends_as_annotated", False),
            )
            print(f"       Done in {time.time()-t0:.1f}s")

        segments = Path(database_dir) / "chimeric_segments.parquet"
        if (getattr(args, "fusion_calls", False) or args.all_post) and segments.exists():
            print("[Post] Calling fusions from chimeric reads...")
            t0 = time.time()
            rescue_bams, rescue_note = _fusion_rescue_bams(args, database_dir)
            if rescue_note:
                print(f"       Fusion rescue (reads the aligner did not split): {'partly, ' if rescue_bams else 'not run, '}{rescue_note}")
            tsv, listed, n_pass = isoforge_core.post_fusions(
                db_handle, str(out_prefix),
                min_reads=args.fusion_min_reads,
                min_reads_intraexonic=args.fusion_min_reads_intraexonic,
                min_mapq=args.fusion_min_mapq,
                readthrough_distance=args.fusion_readthrough_distance,
                min_fraction=args.fusion_min_fraction,
                max_partners=args.fusion_max_partners,
                single_dataset_min_reads=args.fusion_single_dataset_min_reads,
                foldback_min_reads=getattr(args, "fusion_foldback_min_reads", 5),
                foldback_gap=getattr(args, "fusion_foldback_gap", 1000),
                hub_dominance=getattr(args, "fusion_hub_dominance", 20.0),
                hub_dominance_min_reads=getattr(args, "fusion_hub_dominance_min_reads", 10),
                genome_fasta=getattr(args, "genome_fasta", None),
                rescue_bams=rescue_bams,
                pseudogene_partners=getattr(args, "fusion_pseudogene_partners", True),
                unnamed_partners=getattr(args, "fusion_unannotated_partners", True),
                min_reads_unnamed=getattr(args, "fusion_min_reads_unannotated", 5),
                min_reads_unnamed_intraexonic=2 * getattr(args, "fusion_min_reads_unannotated", 5),
            )
            print(f"       {n_pass} fusions pass of {listed} candidates -> {tsv} ({time.time()-t0:.1f}s)")

        if getattr(args, "hypergraph", False):
            print("[Post] Building the read hypergraph and testing splicing coordination...")
            t0 = time.time()
            cpath, n_edges, n_pairs, n_sig = isoforge_core.post_hypergraph(
                db_handle, str(out_prefix),
                whitelist_file=str(whitelist_file) if whitelist_file else None,
                min_reads=getattr(args, "coordination_min_reads", 20.0),
            )
            print(f"       {n_edges} hyperedges; {n_sig} of {n_pairs} pairs of alternative segments coordinated (q < 0.05) -> {cpath} ({time.time()-t0:.1f}s)")

        if args.read_annot or args.all_post:
            print("[Post] Generating read annotations...")
            t0 = time.time()
            isoforge_core.post_read_annot(
                db_handle, str(out_prefix), args.genome,
                whitelist_file=str(whitelist_file) if whitelist_file else None,
                split_unsupported=getattr(args, "split_unsupported_reads", False),
                max_frac_a=args.max_frac_a_filter,
                frac_a_annotated_end_window=args.frac_a_annotated_end_window,
                emit_tsv=getattr(args, "read_annot_format", "both") != "parquet",
            )
            print(f"       Done in {time.time()-t0:.1f}s")

        if whitelist_file and Path(whitelist_file).exists():
            print("[Post] Pruning cluster BED/Parquet to whitelisted transcripts...")
            isoforge_core.post_prune_clusters(db_handle, str(whitelist_file))
    finally:
        # Always free the Rust-allocated PostDB
        isoforge_core.post_free_db(db_handle)


def main():
    args = get_args()

    # ── Import Rust engine ────────────────────────────────────────────────
    try:
        import isoforge_core
    except ImportError:
        print(
            "ERROR: the isoforge_core Rust extension is not installed.\n"
            "  Install IsoForge from its repository: bash install.sh (or pip install .)",
            file=sys.stderr,
        )
        sys.exit(1)
    _engine = isoforge_core.version()
    print(f"[IsoForge] IsoForge annotate v{__version__}" + (f" (engine {_engine})" if _engine != __version__ else ""), flush=True)

    # ── Output paths ──────────────────────────────────────────────────────
    # --o is an output prefix; its directory is created if needed.
    # e.g. --o /results/project → dir=/results, prefix=project
    # e.g. --o project           → dir=., prefix=project
    out_path = Path(args.outprefix)
    # A trailing "_isoforge" is stripped because output names append it.
    stem = out_path.name
    if stem.endswith("_isoforge"):
        stem = stem[: -len("_isoforge")]
        out_path = out_path.parent / stem
    out_dir = out_path.parent
    out_prefix = str(out_path)           # full path prefix for output files
    project_name = out_path.name         # just the prefix name (e.g. project)
    database_dir = out_dir / (project_name + "_isoforge_db")

    # ── Genome FASTA ──────────────────────────────────────────────────────
    if args.genome_fasta:
        if not Path(args.genome_fasta).exists():
            print(f"ERROR: Genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            sys.exit(1)
        ensure_fai(args.genome_fasta)

    # ── Resolve the samples ───────────────────────────────────────────────
    samples = input_samples(args)
    from_reads = [is_read_file(path) for path, _ in samples]
    if any(from_reads) and not all(from_reads):
        print("ERROR: give every sample as a read file or every sample as a BAM, not a mix.", file=sys.stderr)
        sys.exit(1)
    from_reads = all(from_reads)
    if from_reads:
        if not args.genome_fasta:
            print("ERROR: aligning reads needs the genome FASTA (--genome).", file=sys.stderr)
            sys.exit(1)
        if len({name for _, name in samples}) != len(samples):
            print("ERROR: the sample names must differ from each other (--names).", file=sys.stderr)
            sys.exit(1)

    # ── Resolve annotation directory ──────────────────────────────────────
    annot_dir = args.annot_dir

    if annot_dir is None and args.annotation:
        # The index of the annotation GTF, built when there is none yet
        annot_dir = str(index_for_gtf(args.annotation, options.default_build(args), out_dir))

    if annot_dir is None:
        # Auto-discover isoforge_index/
        candidates = []

        # 1. Next to the genome FASTA
        if args.genome_fasta:
            fasta_parent = Path(args.genome_fasta).resolve().parent
            candidates.append(fasta_parent / "isoforge_index")

        # 2. Current working directory
        candidates.append(Path.cwd() / "isoforge_index")

        for cand in candidates:
            if cand.is_dir() and (cand / "edges.parquet").exists():
                annot_dir = str(cand)
                print(f"[IsoForge] Auto-discovered index at: {cand}")
                break

        if annot_dir is None:
            searched = ", ".join(str(c) for c in candidates)
            print(
                "ERROR: no annotation given and no isoforge_index/ found.\n"
                f"  Searched: {searched}\n"
                "  Give the annotation GTF (--annotation) or an index built by `isoforge index` (--db).",
                file=sys.stderr,
            )
            sys.exit(1)

    if not Path(annot_dir).is_dir():
        print(f"ERROR: IsoForge index directory not found: {annot_dir}", file=sys.stderr)
        print("  Give the annotation GTF (--annotation), or build the index with `isoforge index`.", file=sys.stderr)
        sys.exit(1)

    # The genome build name: as given, else the one the index records
    _index_build = (index_metadata(annot_dir) or {}).get("genome_build")
    if not args.genome:
        if not _index_build:
            print("ERROR: the index records no genome build name; give it with --build NAME.", file=sys.stderr)
            sys.exit(1)
        args.genome = _index_build
    elif _index_build and _index_build != args.genome:
        print(
            f"WARNING: Index genome_build '{_index_build}' does not match "
            f"--genome '{args.genome}'. Proceeding, but results may be incorrect.",
            file=sys.stderr,
        )

    # ── Resolve SJ files ──────────────────────────────────────────────────
    sj_files = args.sj_files or []
    if args.sj_dir:
        sj_files.extend(discover_sj_files(args.sj_dir))

    # ── Create output directory ───────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    database_dir.mkdir(parents=True, exist_ok=True)

    # ── Align the reads, or take the BAMs as they are ─────────────────────
    aligned = []
    if from_reads:
        aligned = align_samples(args, samples, out_prefix)
        bam_configs = [(bam, name) for bam, (_, name) in zip(aligned, samples)]
    else:
        bam_configs = list(samples)
    # the run's BAMs, for the fusion rescue
    args.run_bams = [(name, bam) for bam, name in bam_configs]

    # ── Build config dict ─────────────────────────────────────────────────
    config = {
        "repair_microexons": args.repair_microexons,
        "splice_realign_window": args.splice_realign_window,
        "terminal_exon_rescue": args.terminal_exon_rescue,
        "sequence_settled_sites": args.sequence_settled_sites,
        "approx_chain_tolerance": args.approx_chain_tolerance,
        "chain_first": args.chain_first,
        "approx_chain_lock_annotated": not args.approx_chain_no_lock,
        "assign_unique_ism": args.assign_unique_ism,
        "ism_reassign_max_reads": args.ism_reassign_max_reads,
        "ism_intact_end_support": args.ism_intact_end_support,
        "known_ends_not_degraded": args.known_ends_not_degraded,
        "nested_fsm_ambiguous": args.nested_fsm_ambiguous,
        "tss_peak_raw_mode": args.tss_peak_raw_mode,
        "annotated_id_follows_ends": args.annotated_id_follows_ends,
        "exclude_upstream_unclustered_5p": args.exclude_upstream_unclustered_5p,
        "split_capped_tss": args.split_capped_tss,
        "split_tss_min_distance": args.split_tss_min_distance,
        "split_tss_min_reads": args.split_tss_min_reads,
        "split_tss_min_share": args.split_tss_min_share,
        "split_tss_cap_ratio": args.split_tss_cap_ratio,
        "nested_fsm_end_window": args.nested_fsm_end_window,
        "compact_tss": args.compact_tss,
        "reference_end_support": args.reference_end_support,
        "known_singleton_ends": args.known_singleton_ends,
        "threads": args.threads,
        "batch_size": args.batch_size,
        "min_coverage": args.min_coverage,
        "min_length": args.min_read_length,
        "small_rna": args.small_rna,
        "small_rna_types": [t.strip() for t in args.small_rna_types.split(",") if t.strip()],
        "small_rna_min_read_length": args.small_rna_min_read_length,
        "small_rna_max_length": args.small_rna_max_length,
        "coverage_excludes_polya": args.coverage_excludes_polya,
        "junction_consensus_window": args.junction_consensus_window,
        "novel_site_window": args.novel_site_window,
        "splice_scores": args.splice_scores,
        "min_identity": args.min_identity,
        "max_frac_as": args.max_frac_as,
        "max_5p_dist": args.max_5p,
        "max_3p_dist": args.max_3p,
        "max_splice_dist": args.max_splice_dist,
        "sj_correction_window": args.sj_correction_window,
        "max_extended_snap": args.max_extended_snap,
        "create_novel_spliced_genes": args.create_novel_spliced_genes,
        "fsm_terminal_new_id": args.fsm_terminal_new_id,
        "strict_full_length": args.strict_full_length,
        "infer_strand_from_polya": args.infer_strand_from_polya,
        "id_prefix": args.id_prefix,
        "verbosity": args.verbosity,
        "correct_sj": args.correct_sj,
        "correct_microindels": args.correct_microindels,
        "max_indel_size": args.max_indel_size,
        "min_mapq": args.min_mapq,
        "quant_min_mapq": max(0, min(255, args.quant_min_mapq)),
        "detect_fusions": bool(getattr(args, "fusions", True)),
        "fusion_spliced_joins": bool(getattr(args, "fusion_spliced_joins", False)),
        "junction_hmm": bool(getattr(args, "junction_hmm", False)),
        "junction_hmm_window": int(getattr(args, "junction_hmm_window", 5)),
        "junction_hmm_min_support": int(getattr(args, "junction_hmm_min_support", 2)),
        "junction_hmm_model": getattr(args, "junction_hmm_model", None),
        "junction_indel_min": max(0, args.junction_indel_min),
    }

    if (args.validate_novel_junctions or args.repair_microexons) and not args.genome_fasta:
        print("  NOTE: no genome FASTA (--genome): novel junction motifs are not checked and microexons are not repaired "
              "(--validate_novel_junctions, --repair_microexons)")
        args.validate_novel_junctions = False
        args.repair_microexons = False
        config["repair_microexons"] = False
    if args.small_rna is None:
        args.small_rna = bool(args.genome_fasta)
        if not args.small_rna:
            print("  NOTE: small RNAs are not counted without the genome FASTA (--genome)")
    elif args.small_rna and not args.genome_fasta:
        print("ERROR: --small_rna requires the genome FASTA (--genome)", file=sys.stderr)
        sys.exit(1)
    config["small_rna"] = args.small_rna
    if args.min_count_supported != 0 and args.min_count_supported < 5:
        raise ValueError("--minCountSupported must be 0 (disabled) or at least 5")

    # Integrated preprocessing
    if args.genome_fasta:
        if not Path(args.genome_fasta).exists():
            print(f"ERROR: Genome FASTA not found: {args.genome_fasta}", file=sys.stderr)
            sys.exit(1)
        config["genome_fasta"] = args.genome_fasta

    if sj_files:
        config["sj_files"] = sj_files

    # TSS/TES clustering — on with the genome FASTA unless --no_clustering
    if args.genome_fasta and not args.no_clustering:
        config["clustering"] = True
        config["gaussian_std"] = args.gaussian_std
        config["merge_distance"] = args.merge_distance
        config["min_cluster_reads"] = args.min_cluster_reads
        config["tss_peak_max_entropy"] = args.tss_peak_max_entropy
        config["tss_peak_rescue_min_reads"] = args.tss_peak_rescue_min_reads
        config["tss_peak_rescue_enrichment"] = args.tss_peak_rescue_enrichment
        config["end_contradiction_distance"] = args.end_contradiction_distance
        config["annotated_tes_polya_window"] = args.frac_a_annotated_end_window
        config["reference_end_window"] = args.reference_end_window
        config["min_condition_reads"] = args.min_condition_reads
        config["min_replication_rate"] = args.min_replication_rate
        config["filter_internal_priming"] = args.filter_internal_priming
        config["filter_truncated_reads"] = args.filter_truncated_reads
        config["degradation_tolerance"] = args.degradation_tolerance

        # Orthogonal evidence gates
        config["require_polya_tail"] = args.require_polya_tail
        config["write_end_tracks"] = getattr(args, "write_end_tracks", True)

        if args.cage_peaks:
            cage_path = Path(args.cage_peaks)
            if not cage_path.exists():
                print(f"ERROR: CAGE peaks file not found: {args.cage_peaks}", file=sys.stderr)
                sys.exit(1)
            config["cage_peaks_file"] = str(cage_path)
            config["cage_window"] = args.cage_window

    # ── Banner ────────────────────────────────────────────────────────────
    print(f"{'═' * 60}")
    print(f"  IsoForge v{isoforge_core.version()} — 5′ Isoforms · Ends · Counts · Fusions 3′")
    print(f"{'═' * 60}")
    print(f"  Prefix:        {out_prefix}")
    print(f"  Index:         {annot_dir}")
    print(f"  Genome build:  {args.genome}")
    print(f"  Datasets:      {len(bam_configs)}")
    for bam, name in bam_configs:
        print(f"    • {name}: {bam}")
    print(f"  Threads:       {args.threads}")
    print(f"  Min MAPQ:      {args.min_mapq}")
    print(f"  Output dir:    {out_dir}")
    print(f"  Run database:  {database_dir}{' (will be kept)' if args.keep_db else ' (auto-cleaned after post-processing)'}")
    if args.genome_fasta:
        print(f"  Genome FASTA:  {args.genome_fasta}")
        print(f"  SJ correction: {'enabled' if args.correct_sj else 'disabled'}")
        print(f"  Indel correction: {'enabled (max {0}bp)'.format(args.max_indel_size) if args.correct_microindels else 'disabled'}")
    if sj_files:
        print(f"  SJ files:      {len(sj_files)} file(s)")
    if config.get("clustering"):
        print(f"  TSS/TES clustering: enabled")
        print(f"    σ={args.gaussian_std}  merge={args.merge_distance}bp  "
              f"min_reads={args.min_cluster_reads}  min_cond={args.min_condition_reads}  "
              f"repl_rate={args.min_replication_rate}")
        print(f"    filter_IP={'yes' if args.filter_internal_priming else 'no'}  "
              f"filter_truncated={'yes' if args.filter_truncated_reads else 'no'}")
        print(f"    degradation_tol={args.degradation_tolerance}bp  "
              f"cluster_gate=OR (either TSS or TES)")
        # Orthogonal evidence gates
        print(f"  Poly(A) tail filter: {'ENABLED (mandatory)' if args.require_polya_tail else 'DISABLED'}")
        if config.get("cage_peaks_file"):
            print(f"  CAGE peak filter:   ENABLED")
            print(f"    peaks: {args.cage_peaks}")
            print(f"    window: {args.cage_window}bp")
        else:
            print(f"  CAGE peak filter:   disabled (no --cage_peaks provided)")
    else:
        print(f"  TSS/TES clustering: disabled")
    print(f"{'─' * 60}")

    # ── Enable all post-processing if --all ───────────────────────────────
    any_post = args.summarize or args.filter or args.abundance or args.gtf or args.read_annot or args.all_post

    # ── Initialize engine ─────────────────────────────────────────────────
    print("\n[1/3] Initializing IsoForge engine...")
    t_start = time.time()

    try:
        engine = isoforge_core.IsoForgeEngine(
            annot_dir=annot_dir,
            config=config,
            output_dir=str(database_dir),
        )
    except Exception as e:
        print(f"ERROR: Engine initialization failed: {e}", file=sys.stderr)
        sys.exit(1)

    t_init = time.time() - t_start
    print(f"       Engine loaded in {t_init:.1f}s")

    # ── Annotate ──────────────────────────────────────────────────────────
    print(f"\n[2/3] Annotating {len(bam_configs)} dataset(s)...")
    t_annot_start = time.time()

    try:
        result = engine.run(bam_configs)
    except Exception as e:
        print(f"ERROR: Annotation failed: {e}", file=sys.stderr)
        sys.exit(1)

    t_annot = time.time() - t_annot_start
    print(f"       Annotation complete in {t_annot:.1f}s")

    # ── Write run-time counters (id_prefix etc.) into results DB ──────────
    import json as _json
    _run_counters = {"id_prefix": args.id_prefix}
    _run_counters_path = database_dir / "counters.json"
    with open(_run_counters_path, "w") as _f:
        _json.dump(_run_counters, _f, indent=2)
    # the run's BAMs, for `isoforge post` to re-read (fusion reads the aligner did not split)
    with open(database_dir / "bam_paths.csv", "w") as _f:
        for _bam, _name in bam_configs:
            _f.write(f"{_name},{os.path.abspath(_bam)}\n")

    # ── Free engine memory before post-processing ─────────────────────────
    del engine, result
    import gc; gc.collect()

    # ── Post-processing ───────────────────────────────────────────────────
    if any_post:
        print(f"\n[3/3] Post-processing...")
        t_post_start = time.time()

        run_post_processing(annot_dir, str(database_dir), out_prefix, args, project_name)

        t_post = time.time() - t_post_start
        print(f"\n       Post-processing complete in {t_post:.1f}s")
    else:
        print(f"\n[3/3] Post-processing skipped")

    # ── Clean up intermediate database if not keeping ─────────────────────
    if not getattr(args, 'keep_db', False) and any_post and database_dir.exists():
        import shutil
        db_size = sum(f.stat().st_size for f in database_dir.rglob("*") if f.is_file())
        shutil.rmtree(database_dir)
        print(f"\n  Cleaned up intermediate database ({db_size / 1_000_000:.1f} MB)")
        print(f"  (Use --keep-db to retain intermediate Parquet files)")

    # ── Remove the alignments made from reads unless they are kept ────────
    if aligned and not args.write_bam:
        remove_alignments(aligned)
        print(f"\n  Removed the alignments (use --write_bam to keep them)")

    # ── Summary ───────────────────────────────────────────────────────────
    t_total = time.time() - t_start
    hours, remainder = divmod(t_total, 3600)
    minutes, seconds = divmod(remainder, 60)

    print(f"\n{'═' * 60}")
    print(f"  IsoForge Complete")
    print(f"{'═' * 60}")
    print(f"  Total runtime:  {int(hours)}h {int(minutes)}m {seconds:.1f}s")
    print(f"  Output prefix:  {out_prefix}")

    # List output files
    print(f"\n  Output files:")
    for f in sorted(out_dir.glob(f"{project_name}*")):
        if f.is_file():
            size = f.stat().st_size
            if size > 1_000_000:
                print(f"    {f.name:<45} {size/1_000_000:.1f} MB")
            else:
                print(f"    {f.name:<45} {size/1_000:.1f} KB")
    if database_dir.exists():
        print(f"\n  Intermediate database: {database_dir}/")
        for f in sorted(database_dir.glob("*")):
            if f.is_file():
                size = f.stat().st_size
                if size > 1_000_000:
                    print(f"    {f.name:<45} {size/1_000_000:.1f} MB")
                else:
                    print(f"    {f.name:<45} {size/1_000:.1f} KB")

    # ── Output manifest (ISOFORGE_OUTPUTS_MANIFEST) ───────────────────────
    if os.environ.get("ISOFORGE_OUTPUTS_MANIFEST"):
        manifest_path = out_dir / f"{project_name}_isoforge_outputs.json"
        try:
            import json as _json

            def _entry(path):
                return {"name": path.name, "path": str(path.resolve()),
                        "bytes": path.stat().st_size}

            manifest = {
                "isoforge_version": __version__,
                "project": project_name,
                "prefix": str(out_prefix),
                "output_dir": str(out_dir.resolve()),
                "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "files": [_entry(f) for f in sorted(out_dir.glob(f"{project_name}*"))
                          if f.is_file() and f.name != manifest_path.name],
                "directories": {
                    d.name: sorted(x.name for x in d.iterdir())
                    for d in sorted(out_dir.glob(f"{project_name}*")) if d.is_dir()
                },
                "database_kept": database_dir.exists(),
                "database_dir": str(database_dir.resolve()) if database_dir.exists() else None,
            }
            manifest_path.write_text(_json.dumps(manifest, indent=2) + "\n")
            print(f"\n  Manifest: {manifest_path.name}")
        except Exception as exc:  # an error writing the manifest does not stop the run
            print(f"\n  NOTE: could not write output manifest: {exc}")

    print(f"\n  ✓ Done.")


if __name__ == "__main__":
    main()
