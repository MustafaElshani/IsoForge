#!/usr/bin/env python3
"""isoforge post — filter, count and write the outputs again from a kept run.

Opens the run database kept by `isoforge annotate --keep-db` and repeats the steps that follow annotation.

    isoforge post --database results/project_isoforge_db --annotation gencode.v49.gtf \\
                  --genome hg38.fa --minCount 10 --o results/project_min10

It takes the filtering options of `isoforge annotate`.
"""

import sys
from pathlib import Path

from isoforge import options
from isoforge.genome_fasta import ensure_fai
from isoforge.isoforge_init import index_for_gtf, index_metadata

__version__ = "1.0.0"


POST_EPILOG = """\
Example:
  isoforge post --database results/project_isoforge_db --annotation gencode.v49.gtf \\
      --genome hg38.fa --minCount 10 --o results/project_min10

Give --genome and --splice_scores as in the run. `isoforge post --help-all` lists every option.
"""


def get_args(argv=None):
    given = sys.argv[1:] if argv is None else argv
    p = options.parser_for(
        "post",
        usage="isoforge post --database RUN_DB (--annotation GTF | --db DIR) --genome FASTA --o PREFIX [options]",
        description="Filter, count and write the outputs again from a run database kept with\n"
                    "`isoforge annotate --keep-db`, without reading the alignments again.",
        epilog=POST_EPILOG,
        show_all="--help-all" in given,
        version=f"isoforge post {__version__}",
    )
    args = p.parse_args(argv)

    options.take_genome_fasta(args)
    options.name_outputs(args, ("summarize", "filter", "abundance", "gtf", "fusion_calls", "hypergraph", "read_annot",
                                "all_post"))
    if args.genome_fasta:
        if not Path(args.genome_fasta).exists():
            p.error(f"genome FASTA not found: {args.genome_fasta}")
        ensure_fai(args.genome_fasta)
    if args.validate_novel_junctions and not args.genome_fasta:
        print("[Post] no genome FASTA (--genome): novel junction motifs are not checked (--validate_novel_junctions)")
        args.validate_novel_junctions = False
    if not args.evidence_tss_rescue:
        args.ism_evidence_rescue = False
    if args.min_count_supported != 0 and args.min_count_supported < 5:
        p.error("--minCountSupported must be 0 (disabled) or at least 5")

    if not args.annot_dir:
        if not args.annotation:
            p.error("give the annotation of the run: --annotation GTF, or its index with --db")
        args.annot_dir = str(index_for_gtf(args.annotation, options.default_build(args), Path(args.outprefix).parent))
    for label, path in [("index", args.annot_dir), ("database", args.database_dir)]:
        if not Path(path).is_dir():
            p.error(f"{label} directory not found: {path}")
    if not args.genome:
        # the genome build name the index records
        args.genome = (index_metadata(args.annot_dir) or {}).get("genome_build")
        if not args.genome:
            p.error("the index records no genome build name; give it with --build NAME")

    return args


def main():
    args = get_args()

    # Imported here so `--help` works without a built extension.
    from isoforge.isoforge_main import run_post_processing

    out_prefix = Path(args.outprefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    print(f"IsoForge post v{__version__}")
    print(f"  index    : {args.annot_dir}")
    print(f"  database : {args.database_dir}")
    print(f"  output   : {out_prefix}")

    run_post_processing(Path(args.annot_dir), Path(args.database_dir),
                        out_prefix, args, out_prefix.name)
    print("\nDone.")


if __name__ == "__main__":
    sys.exit(main())
