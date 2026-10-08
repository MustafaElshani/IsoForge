"""The options of `isoforge annotate` and `isoforge post`, and the help layout of every command.

One table serves both commands. `--help` lists the options of everyday use, `--help-all` every option.
"""
import argparse
import os
import shutil

ADMISSION_CLASSES = "NIC:0.5:0.1,NNC:never:0.7"
SIBLING_EDIT_RATIOS = "1,0.5,0.1,0.05"

A, P, AP = ("annotate",), ("post",), ("annotate", "post")


class Switch(argparse.BooleanOptionalAction):
    """An on/off option (--name / --no-name) whose help is shown as written."""

    def __init__(self, option_strings, dest, help=None, **kwargs):
        super().__init__(option_strings, dest, help=help, **kwargs)
        self.help = help


def opt(*flags, cmds=AP, main=(), **settings):
    """An option: its flags, the commands that take it, the commands whose short help lists it, and its settings."""
    return flags, cmds, main, settings


FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def take_genome_fasta(args):
    """Moves a FASTA named by --genome to args.genome_fasta; args.genome keeps the build name given, or None."""
    if args.genome and args.genome.lower().endswith(FASTA_SUFFIXES):
        if not args.genome_fasta:
            args.genome_fasta = args.genome
        args.genome = None
    if args.build:
        args.genome = args.build


def default_build(args):
    """The genome build name of a new index: the name given, else the FASTA's file name up to its first dot."""
    if args.genome:
        return args.genome
    return os.path.basename(args.genome_fasta).split(".")[0] if args.genome_fasta else "genome"


def name_outputs(args, steps):
    """Asks for every output when none of the single outputs is named."""
    if not any(getattr(args, step) for step in steps):
        args.all_post = True


# Each group is a section of the help.
GROUPS = [
    ("Input", [
        opt("--reads", dest="reads", nargs='+', metavar="FASTQ", cmds=A, main=A,
            help="read files (FASTQ or FASTA, may be gzipped), one sample each"),
        opt("--bam", dest="bam", nargs='+', metavar="BAM", cmds=A, main=A,
            help="sorted BAM files, one sample each"),
        opt("--names", dest="names", nargs='+', metavar="NAME", cmds=A, main=A,
            help="sample names for --reads or --bam (default: the file names)"),
        opt("--config", "--f", dest="config_file", metavar="CSV", cmds=A, main=A,
            help="sample sheet, one sample per line: name,description,platform,file"),
        opt("--dataset_name", dest="dataset_name", metavar="NAME", cmds=A,
            help="sample name of a single --bam; as --names"),
        opt("--database", dest="database_dir", required=True, metavar="DIR", cmds=P, main=P,
            help="run database kept by `isoforge annotate --keep-db`"),
    ]),
    ("Reference", [
        opt("--genome", dest="genome", metavar="FASTA", cmds=A, main=A,
            help="genome FASTA"),
        opt("--genome", dest="genome", metavar="FASTA", cmds=P, main=P,
            help="genome FASTA, as given to the run"),
        opt("--annotation", dest="annotation", metavar="GTF", cmds=A, main=A,
            help="annotation GTF; its index is built once and used again"),
        opt("--annotation", dest="annotation", metavar="GTF", cmds=P, main=P,
            help="annotation GTF, as given to the run"),
        opt("--db", "--annot", dest="annot_dir", metavar="DIR", cmds=AP, main=AP,
            help="IsoForge index, in place of --annotation"),
        opt("--genome_fasta", dest="genome_fasta", metavar="FASTA", cmds=AP,
            help="genome FASTA, when --genome gives the build name of the index"),
        opt("--build", dest="build", metavar="NAME", cmds=AP,
            help="genome build name of the index (default: the one the index records)"),
        opt("--splice_scores", dest="splice_scores", metavar="TSV", cmds=A, main=A,
            help="splice-site scores (minisplice table) for new splice sites"),
        opt("--splice_scores", dest="splice_scores", metavar="TSV", cmds=P, main=P,
            help="splice-site scores (minisplice table), as given to the run"),
        opt("--sj_dir", dest="sj_dir", metavar="DIR", cmds=AP, main=A,
            help="directory of STAR SJ.out.tab files: short-read junction evidence"),
        opt("--sj_files", dest="sj_files", nargs='+', metavar="FILE", cmds=AP,
            help="STAR SJ.out.tab files, in place of --sj_dir"),
        opt("--cage_peaks", dest="cage_peaks", metavar="BED", cmds=A, main=A,
            help="CAGE peaks (BED); TSS clusters away from a peak are removed"),
        opt("--cage_window", dest="cage_window", type=int, default=50, metavar="N", cmds=A,
            help="bp between a TSS cluster and a CAGE peak (default: 50)"),
    ]),
    ("Output", [
        opt("--o", dest="outprefix", required=True, metavar="PREFIX", cmds=AP, main=AP,
            help="prefix of the output files; may include a directory"),
        opt("--all", dest="all_post", action="store_true", cmds=AP,
            help="write every output (the default when no single output is named)"),
        opt("--summarize", dest="summarize", action="store_true", cmds=AP,
            help="write the summary table"),
        opt("--filter", dest="filter", action="store_true", cmds=AP,
            help="write the list of transcripts that pass the filters"),
        opt("--abundance", dest="abundance", action="store_true", cmds=AP,
            help="write the abundance table"),
        opt("--gtf", dest="gtf", action="store_true", cmds=AP,
            help="write the GTF of the transcripts that pass"),
        opt("--read_annot", dest="read_annot", action="store_true", cmds=AP,
            help="write the read table"),
        opt("--fusions", dest="fusion_calls", action="store_true", cmds=P,
            help="write the fusion calls"),
        opt("--hypergraph", dest="hypergraph", action="store_true", cmds=AP, main=AP,
            help="write the read hypergraph and the splicing-coordination tests"),
        opt("--read_annot_format", dest="read_annot_format", choices=('both', 'parquet', 'tsv'), default='both',
            cmds=AP, main=AP,
            help="read table as Parquet, TSV or both (default: both)"),
        opt("--unfiltered", dest="unfiltered", action="store_true", cmds=A, main=A,
            help="also write the abundance table before filtering"),
        opt("--keep-db", dest="keep_db", action="store_true", cmds=A, main=A,
            help="keep the run database, for `isoforge post`"),
        opt("--end_tracks", dest="write_end_tracks", action=Switch, default=True, cmds=A, main=A,
            help="write read-end coverage tracks per dataset (default: on)"),
        opt("--split_unsupported_reads", dest="split_unsupported_reads", action=Switch, default=True, cmds=AP,
            help="write reads that are not whole reads of their model to a table of their own (default: on)"),
        opt("--coordination_min_reads", dest="coordination_min_reads", type=float, default=20.0, metavar="N",
            cmds=AP,
            help="reads spanning both segments that a coordination test needs (default: 20)"),
        opt("--idprefix", dest="id_prefix", default='ISOFOR', metavar="PREFIX", cmds=A,
            help="prefix of novel gene and transcript IDs (default: ISOFOR)"),
    ]),
    ("Alignment of reads", [
        opt("--platform", dest="platform", choices=('ont', 'hifi'), default='ont', cmds=A, main=A,
            help="ont: minimap2 -x splice; hifi: -x splice:hq (default: ont)"),
        opt("--oriented", dest="oriented", action="store_true", cmds=A, main=A,
            help="reads are in transcript orientation, e.g. direct RNA (minimap2 -uf)"),
        opt("--copy_comment", dest="copy_comment", action="store_true", cmds=A, main=A,
            help="keep FASTQ comments as tags, e.g. poly(A) tail lengths (minimap2 -y)"),
        opt("--write_bam", dest="write_bam", action="store_true", cmds=A, main=A,
            help="keep the alignments in <prefix>_alignments/"),
        opt("-k", dest="k", type=int, metavar="N", cmds=A,
            help="minimizer length (default: the preset's)"),
        opt("--junc_bed", dest="junc_bed", metavar="BED", cmds=A,
            help="annotated junctions for the aligner (minimap2 --junc-bed)"),
        opt("--splice_score_scale", dest="splice_score_scale", type=float, default=0.7, metavar="X", cmds=A,
            help="weight of --splice_scores in the alignment, which uses them without --junc_bed (default: 0.7)"),
        opt("--secondary", dest="secondary", action="store_true", cmds=A,
            help="keep secondary alignments"),
        opt("--mm_index", dest="mm_index", metavar="MMI", cmds=A,
            help="minimap2 index: read when the file exists, else built and saved there"),
        opt("--no_sort", dest="no_sort", action="store_true", cmds=A,
            help="leave the alignments unsorted; annotation expects sorted alignments"),
        opt("--sort_mem_gb", dest="sort_mem_gb", type=float, default=4.0, metavar="X", cmds=A,
            help="memory for sorting, in GB (default: 4)"),
        opt("--align_threads", dest="align_threads", type=int, metavar="N", cmds=A,
            help="threads of the alignment (default: --threads)"),
    ]),
    ("Reads", [
        opt("--mapq", dest="min_mapq", type=int, default=16, metavar="N", cmds=A, main=A,
            help="lowest mapping quality of reads that build models (default: 16)"),
        opt("--quant_min_mapq", dest="quant_min_mapq", type=int, default=1, metavar="N", cmds=A,
            help="reads from this mapping quality up to --mapq are counted but build no model (default: 1)"),
        opt("--min_read_length", dest="min_read_length", type=int, default=200, metavar="N", cmds=A, main=A,
            help="shortest read kept, in bases (default: 200)"),
        opt("--min_coverage", "--cov", dest="min_coverage", type=float, default=0.8, metavar="X", cmds=A,
            help="lowest aligned fraction of a read (default: 0.8)"),
        opt("--min_identity", "--identity", dest="min_identity", type=float, default=0.0, metavar="X", cmds=A,
            help="lowest alignment identity (default: 0)"),
        opt("--coverage_excludes_polya", dest="coverage_excludes_polya", action=Switch, default=True, cmds=A,
            help="leave a poly(A) soft clip out of the aligned fraction (default: on)"),
        opt("--no_infer_strand_from_polya", dest="infer_strand_from_polya", action="store_false", cmds=A,
            help="take each read's strand from its alignment alone, not from its poly(A) tail or splice motifs"),
        opt("--no_polya_filter", dest="require_polya_tail", action="store_false", cmds=A, main=A,
            help="accept 3' ends without a poly(A) tail (reads without tail data)"),
        opt("--max_frac_as", dest="max_frac_as", type=float, default=0.5, metavar="X", cmds=A,
            help="A fraction of the 20 nt after a 3' end above which a read without a poly(A) signal places no TES "
                 "(default: 0.5)"),
        opt("--maxFracA", dest="max_frac_a_filter", type=float, default=0.35, metavar="X", cmds=AP,
            help="A fraction above which such a read is internally primed and not counted (default: 0.35)"),
        opt("--frac_a_annotated_end_window", dest="frac_a_annotated_end_window", type=int, default=30, metavar="N",
            cmds=AP,
            help="a read ending within this many bp of an annotated TES is never internally primed (default: 30)"),
    ]),
    ("Splice junctions", [
        opt("--no_correct_sj", dest="correct_sj", action="store_false", cmds=A,
            help="do not correct splice junctions"),
        opt("--max_splice_dist", dest="max_splice_dist", type=int, default=6, metavar="N", cmds=A,
            help="bp a splice site may move to a canonical motif or to an existing site (default: 6)"),
        opt("--sj_correction_window", dest="sj_correction_window", type=int, metavar="N", cmds=A,
            help="search radius of the motif correction, in bp (default: --max_splice_dist; 0 disables)"),
        opt("--no_correct_indels", dest="correct_microindels", action="store_false", cmds=A,
            help="do not absorb small insertions and deletions beside junctions"),
        opt("--max_indel_size", dest="max_indel_size", type=int, default=5, metavar="N", cmds=A,
            help="largest insertion or deletion absorbed, in bp (default: 5)"),
        opt("--junction_indel_min", dest="junction_indel_min", type=int, default=3, metavar="N", cmds=A,
            help="an intron with an indel of this size or more beside it is placed where the read needs none "
                 "(default: 3; 0 disables)"),
        opt("--junction_hmm", dest="junction_hmm", action=Switch, default=True, cmds=A,
            help="settle junctions whose reads disagree by a few bases with a pair-HMM over the reads (default: "
                 "on)"),
        opt("--junction_hmm_window", dest="junction_hmm_window", type=int, default=5, metavar="N", cmds=A,
            help="bp within which introns form one junction group (default: 5)"),
        opt("--junction_hmm_min_support", dest="junction_hmm_min_support", type=int, default=2, metavar="N", cmds=A,
            help="reads an unannotated intron needs to be a candidate of its group (default: 2)"),
        opt("--junction_hmm_model", dest="junction_hmm_model", metavar="JSON", cmds=A,
            help="pair-HMM parameters to use in place of those trained on the run"),
        opt("--junction_consensus_window", dest="junction_consensus_window", type=int, default=6, metavar="N",
            cmds=A,
            help="a new splice site is placed at the most used position within this many bp (default: 6; 0 "
                 "disables)"),
        opt("--novel_site_window", dest="novel_site_window", type=int, default=20, metavar="N", cmds=A,
            help="with --splice_scores, a new splice site moves to the strongest site within this many bp "
                 "(default: 20; 0 disables)"),
        opt("--max_extended_snap", dest="max_extended_snap", type=int, default=25, metavar="N", cmds=A,
            help="a new splice site within this many bp of an annotated one moves onto it unless short reads "
                 "support it (default: 25; 0 disables)"),
        opt("--splice_realign_window", dest="splice_realign_window", type=int, default=30, metavar="N", cmds=A,
            help="an unannotated intron sharing a site with an annotated intron within this many bp is tested "
                 "against it on the read's bases (default: 30; 0 disables)"),
        opt("--sequence_settled_sites", dest="sequence_settled_sites", action=Switch, default=False, cmds=A,
            help="keep a new splice site that the read's bases favour by three edits or more (default: off)"),
        opt("--terminal_exon_rescue", dest="terminal_exon_rescue", action=Switch, default=True, cmds=A,
            help="recover a short first or last exon from the soft clip (default: on)"),
        opt("--repair_microexons", dest="repair_microexons", action=Switch, default=True, cmds=A,
            help="restore annotated microexons the aligner skipped, where the read's bases match (default: on)"),
    ]),
    ("Matching reads to transcripts", [
        opt("--approx_chain_tolerance", dest="approx_chain_tolerance", type=int, default=50, metavar="N", cmds=A,
            help="a read whose splice sites all lie within this many bp of an annotated transcript's is assigned "
                 "to it (default: 50; 0 disables)"),
        opt("--chain_first", dest="chain_first", action=Switch, default=True, cmds=A,
            help="try that whole-chain assignment before making a new combination (default: on)"),
        opt("--approx_chain_no_lock", dest="approx_chain_no_lock", action="store_true", cmds=A,
            help="let whole-chain assignment move a site that lies exactly on an annotated site"),
        opt("--assign_unique_ism", dest="assign_unique_ism",
            choices=('off', 'always', 'below_novel_floor', 'no_fsm_parent'), default='below_novel_floor', cmds=A,
            help="give a truncated read with one compatible annotated transcript to that transcript (default: "
                 "below_novel_floor)"),
        opt("--ism_reassign_max_reads", dest="ism_reassign_max_reads", type=int, default=9, metavar="N", cmds=A,
            help="largest truncated model dissolved under below_novel_floor, in reads (default: 9)"),
        opt("--nested_fsm_ambiguous", dest="nested_fsm_ambiguous", action="store_true", cmds=A,
            help="a read of a short isoform nested in a longer one does not establish the short isoform alone"),
        opt("--nested_fsm_end_window", dest="nested_fsm_end_window", type=int, default=100, metavar="N", cmds=A,
            help="bp from the short isoform's annotated end beyond which such a read is ambiguous (default: 100)"),
        opt("--max_5p", "--5p", dest="max_5p", type=int, default=500, metavar="N", cmds=A,
            help="largest distance between a read's 5' end and an annotated start it matches, in bp (default: 500)"),
        opt("--max_3p", "--3p", dest="max_3p", type=int, default=200, metavar="N", cmds=A,
            help="largest distance between a read's 3' end and an annotated end it matches, in bp (default: 200)"),
        opt("--fsm_terminal_new_id", dest="fsm_terminal_new_id", action="store_true", cmds=A,
            help="give a full match with ends beyond those distances a new ID"),
        opt("--create_novel_spliced_genes", "--nsg", dest="create_novel_spliced_genes", action="store_true", cmds=A,
            help="create novel genes for intergenic spliced transcripts"),
    ]),
    ("Transcript ends", [
        opt("--no_clustering", dest="no_clustering", action="store_true", cmds=A,
            help="do not call TSS and TES clusters"),
        opt("--gaussian_std", dest="gaussian_std", type=float, default=1.0, metavar="X", cmds=A,
            help="standard deviation of the Gaussian that smooths read-end counts (default: 1)"),
        opt("--merge_distance", dest="merge_distance", type=int, default=25, metavar="N", cmds=A,
            help="clusters with peaks within this many bp are merged (default: 25)"),
        opt("--min_cluster_reads", dest="min_cluster_reads", type=int, default=3, metavar="N", cmds=A,
            help="reads a cluster needs (default: 3)"),
        opt("--min_replication_rate", dest="min_replication_rate", type=float, default=0.5, metavar="X", cmds=A,
            help="fraction of the datasets that must hold a cluster, with two or more datasets (default: 0.5; 0 "
                 "disables)"),
        opt("--min_condition_reads", dest="min_condition_reads", type=int, default=5, metavar="N", cmds=A,
            help="a cluster below --min_replication_rate is kept with this many reads in all (default: 5)"),
        opt("--no_filter_internal_priming", dest="filter_internal_priming", action="store_false", cmds=A,
            help="let internally primed reads place TES clusters"),
        opt("--no_filter_truncated_reads", dest="filter_truncated_reads", action="store_false", cmds=A,
            help="let truncated reads place TSS and TES clusters"),
        opt("--degradation_tolerance", dest="degradation_tolerance", type=int, default=20, metavar="N", cmds=A,
            help="bp from its model's most used end beyond which a read end is degraded (default: 20)"),
        opt("--tss_peak_raw_mode", dest="tss_peak_raw_mode", action=Switch, default=True, cmds=A,
            help="judge a TSS cluster at its most used base (default: on)"),
        opt("--tss_peak_max_entropy", dest="tss_peak_max_entropy", type=float, default=0.85, metavar="X", cmds=A,
            help="largest positional entropy of a TSS peak (default: 0.85; 1.01 disables)"),
        opt("--tss_peak_rescue_min_reads", dest="tss_peak_rescue_min_reads", type=int, default=10, metavar="N",
            cmds=A,
            help="reads a gene's strongest rejected TSS peak needs to be kept (default: 10; 0 disables)"),
        opt("--tss_peak_rescue_enrichment", dest="tss_peak_rescue_enrichment", type=float, default=5.0,
            metavar="X", cmds=A,
            help="read density such a peak needs over the 200 bp downstream of it (default: 5)"),
        opt("--compact_tss", dest="compact_tss", action="store_true", cmds=A,
            help="seed TSS clusters from compact starts of several reads"),
        opt("--reference_end_support", dest="reference_end_support", action=Switch, default=True, cmds=A,
            help="reads carrying a whole annotated chain confirm its annotated ends (default: on)"),
        opt("--known_singleton_ends", dest="known_singleton_ends", action=Switch, default=True, cmds=A,
            help="one such read is enough (default: on)"),
        opt("--ism_intact_end_support", dest="ism_intact_end_support", action=Switch, default=True, cmds=A,
            help="a truncated read given to its annotated parent supports the end it kept (default: on)"),
        opt("--reference_end_window", dest="reference_end_window", type=int, default=100, metavar="N", cmds=A,
            help="bp between a read's end and an annotated end for the read to confirm it (default: 100)"),
        opt("--known_ends_not_degraded", dest="known_ends_not_degraded", action=Switch, default=True, cmds=A,
            help="a read ending within --reference_end_window of an annotated end is not degraded there (default: "
                 "on)"),
        opt("--end_contradiction_distance", dest="end_contradiction_distance", type=int, default=100, metavar="N",
            cmds=A,
            help="a read group ending more than this many bp outside an annotated transcript does not take its ID "
                 "(default: 100)"),
        opt("--annotated_id_follows_ends", dest="annotated_id_follows_ends", action=Switch, default=True, cmds=A,
            help="the annotated ID goes to the read group whose ends match the annotated ends (default: on)"),
        opt("--exclude_upstream_unclustered_5p", dest="exclude_upstream_unclustered_5p", action=Switch,
            default=True, cmds=A,
            help="a read starting upstream of its model's TSS without a TSS cluster is not counted (default: on)"),
        opt("--split_capped_tss", dest="split_capped_tss", action=Switch, default=True, cmds=A,
            help="keep a capped start near a model's main start as a model of its own (default: on)"),
        opt("--split_tss_min_distance", dest="split_tss_min_distance", type=int, default=20, metavar="N", cmds=A,
            help="bp between such a start and the main one (default: 20)"),
        opt("--split_tss_min_reads", dest="split_tss_min_reads", type=int, default=50, metavar="N", cmds=A,
            help="reads such a start needs (default: 50)"),
        opt("--split_tss_min_share", dest="split_tss_min_share", type=float, default=0.1, metavar="X", cmds=A,
            help="share of the model's reads such a start needs (default: 0.1)"),
        opt("--split_tss_cap_ratio", dest="split_tss_cap_ratio", type=float, default=0.7, metavar="X", cmds=A,
            help="cap-signature rate such a start needs, relative to the main start's (default: 0.7)"),
        opt("--strict_full_length", dest="strict_full_length", action="store_true", cmds=A,
            help="count only the reads that pass end validation"),
        opt("--known_end_snap", dest="known_end_snap", type=int, default=25, metavar="N", cmds=AP,
            help="write the annotated end when the read-derived end lies within this many bp of it (default: 25; 0 "
                 "disables)"),
        opt("--known_end_max_extension", dest="known_end_max_extension", type=int, default=25, metavar="N", cmds=AP,
            help="a read-derived end more than this many bp outside the annotated span needs "
                 "--known_end_extension_min_fraction of the reads (default: 25; 0 disables)"),
        opt("--known_end_extension_min_fraction", dest="known_end_extension_min_fraction", type=float, default=0.5,
            metavar="X", cmds=AP,
            help="share of the end-assigned reads such an end needs (default: 0.5)"),
        opt("--known_ends_as_annotated", dest="known_ends_as_annotated", action=Switch, default=False, cmds=AP,
            main=AP,
            help="write annotated transcripts at their annotated ends (default: off)"),
    ]),
    ("Annotated transcripts", [
        opt("--minCountKnown", dest="min_count_known", type=int, default=2, metavar="N", cmds=AP, main=AP,
            help="reads per dataset an annotated transcript needs (default: 2)"),
        opt("--known_phased_reads", dest="known_phased_reads", type=int, default=1, metavar="N", cmds=AP,
            help="or this many reads with both ends validated, in any datasets (default: 1; 0 disables)"),
        opt("--known_end_evidence", dest="known_end_evidence", type=int, default=2, metavar="N", cmds=AP,
            help="reads that must anchor each annotated end, where reads carry poly(A) tail lengths (default: 2; 0 "
                 "disables)"),
        opt("--known_end_evidence_window", dest="known_end_evidence_window", type=int, default=50, metavar="N",
            cmds=AP,
            help="window of --known_end_evidence, in bp (default: 50)"),
        opt("--filter_known", dest="filter_known", action="store_true", cmds=AP,
            help="hold annotated transcripts to --minCount and --minDatasets, like novel ones"),
    ]),
    ("Novel transcripts", [
        opt("--minCount", dest="min_count", type=int, default=5, metavar="N", cmds=AP, main=AP,
            help="reads per dataset a novel transcript needs (default: 5)"),
        opt("--minDatasets", dest="min_datasets", type=int, default=1, metavar="N", cmds=AP, main=AP,
            help="datasets a novel transcript must pass in (default: 1; 0 = all)"),
        opt("--pool_datasets", dest="pool_datasets", action=Switch, default=False, cmds=AP,
            help="count a novel transcript's reads over all datasets; for replicates of one sample (default: off)"),
        opt("--minCountNIC", dest="min_count_nic", type=int, default=5, metavar="N", cmds=AP,
            help="read floor of new combinations of annotated junctions (default: 5; 0 means --minCount)"),
        opt("--minCountEndVariant", dest="min_count_end_variant", type=int, default=3, metavar="N", cmds=AP,
            help="read floor of an end variant of an annotated transcript (default: 3; 0 disables)"),
        opt("--end_variant_min_datasets", dest="end_variant_min_datasets", type=int, default=0, metavar="N",
            cmds=AP,
            help="datasets in which --minCountEndVariant must hold (default: 0, two thirds of them)"),
        opt("--minCountSupported", dest="min_count_supported", type=int, default=0, metavar="N", cmds=AP,
            help="read floor of compact, end-supported combinations of annotated introns, 5 or more (default: 0, "
                 "off)"),
        opt("--weak_source_tier", dest="weak_source_tier", action=Switch, default=True, cmds=AP,
            help="a new combination passes on 3 reads where the nearest reported model of its gene has at most 5 "
                 "(default: on)"),
        opt("--min_novel_locus_fraction", dest="min_novel_locus_fraction", type=float, default=0.005, metavar="X",
            cmds=AP,
            help="reads a novel transcript needs as a fraction of its gene's top model (default: 0.005; 0 "
                 "disables)"),
        opt("--novel_exact_support", dest="novel_exact_support", type=int, default=3, metavar="N", cmds=AP,
            help="reads per dataset whose splice sites match the model's exactly (default: 3; 0 disables)"),
        opt("--novel_exact_tol", dest="novel_exact_tol", type=int, default=0, metavar="N", cmds=AP,
            help="tolerance of --novel_exact_support, in bp (default: 0)"),
        opt("--validate_novel_junctions", dest="validate_novel_junctions", action=Switch, default=True, cmds=AP,
            help="unannotated introns need GT-AG or GC-AG ends and, with --splice_scores, scored sites (default: "
                 "on)"),
        opt("--junction_prune_radius", dest="junction_prune_radius", type=int, default=50, metavar="N", cmds=AP,
            help="remove a novel model whose intron is a minor variant of another within this many bp (default: "
                 "50; 0 disables)"),
        opt("--junction_prune_ratio", dest="junction_prune_ratio", type=float, default=0.7, metavar="X", cmds=AP,
            help="read ratio below which --junction_prune_radius removes (default: 0.7)"),
        opt("--sibling_edits", dest="sibling_edits", action=Switch, default=True, cmds=AP,
            help="remove a novel model one alignment edit away from a reported model with more reads (default: on)"),
        opt("--sibling_edit_ratios", dest="sibling_edit_ratios", default=SIBLING_EDIT_RATIOS, metavar="A,B,C,D",
            cmds=AP,
            help="read ratios of the four edits: lost terminal exon, skipped microexon, moved junction, extra "
                 "intron (default: 1,0.5,0.1,0.05)"),
        opt("--sub_path_min_fraction", dest="sub_path_min_fraction", type=float, default=0.3, metavar="X", cmds=AP,
            help="share of a longer model's reads that a model contained in it needs (default: 0.3; 0 disables)"),
        opt("--sub_path_tail_evidence", dest="sub_path_tail_evidence", action=Switch, default=True, cmds=AP,
            help="a contained model whose reads carry poly(A) tails at a 3' end of its own passes (default: on)"),
        opt("--minTSSSupport", dest="min_tss_support", type=float, default=0.7, metavar="X", cmds=AP,
            help="share of a novel model's reads that must start in a 5' UTR or intergenic sequence (default: 0.7; "
                 "0 disables)"),
        opt("--tss_upstream_rescue", dest="tss_upstream_rescue", action=Switch, default=True, cmds=AP,
            help="a spliced model with a first exon of its own passes --minTSSSupport (default: on)"),
        opt("--cap_tss_min_fraction", dest="cap_tss_min_fraction", type=float, default=0.15, metavar="X", cmds=AP,
            help="or this share of its clustered starts carries the cap signature (default: 0.15; 0 disables)"),
        opt("--cap_tss_min_reads", dest="cap_tss_min_reads", type=int, default=3, metavar="N", cmds=AP,
            help="capped reads --cap_tss_min_fraction needs (default: 3)"),
        opt("--evidence_tss_rescue", dest="evidence_tss_rescue", action=Switch, default=True, cmds=AP,
            help="compact, end-supported starts outside annotated UTRs pass (default: on)"),
        opt("--ism_evidence_rescue", dest="ism_evidence_rescue", action=Switch, default=True, cmds=AP,
            help="truncated models with a compact start pass the start rule (default: on)"),
        opt("--ism_rescue_min_reads", dest="ism_rescue_min_reads", type=int, default=5, metavar="N", cmds=AP,
            help="reads such a start needs (default: 5)"),
        opt("--ism_rescue_min_distance", dest="ism_rescue_min_distance", type=int, default=0, metavar="N", cmds=AP,
            help="bp such a start must lie from the start of every annotated parent (default: 0)"),
        opt("--require_phased_novel_ends", dest="require_phased_novel_ends", action=Switch, default=True, cmds=AP,
            help="novel unspliced models need 3 reads that join a start and an end cluster (default: on)"),
        opt("--phased_ends_spliced", dest="phased_ends_spliced", action="store_true", cmds=AP,
            help="hold spliced novel models to that rule too"),
        opt("--start_model", dest="start_model", action="store_true", cmds=AP,
            help="remove novel models that start at the acceptor of a longer chain's junction (experimental)"),
        opt("--start_min_llr", dest="start_min_llr", type=float, default=20.0, metavar="X", cmds=AP,
            help="start evidence that exempts a model under --start_model (default: 20)"),
        opt("--allowGenomic", dest="allow_genomic", action=Switch, default=True, cmds=AP,
            help="keep models of the Genomic class (default: on)"),
    ]),
    ("Learned admission", [
        opt("--admission_model", dest="admission_model", metavar="JSON", cmds=AP,
            help="tree ensemble to use (default: the model shipped with IsoForge)"),
        opt("--admission_add", dest="admission_add", type=float, default=0.6, metavar="X", cmds=AP,
            help="probability from which a rejected candidate is added (default: 0.6)"),
        opt("--admission_drop", dest="admission_drop", type=float, default=0.5, metavar="X", cmds=AP,
            help="probability below which a reported model is removed (default: 0.5)"),
        opt("--admission_classes", dest="admission_classes", default=ADMISSION_CLASSES, metavar="LIST", cmds=AP,
            help="thresholds of single classes as CLASS:ADD:VETO, comma-separated; 'none' for none (default: "
                 "NIC:0.5:0.1,NNC:never:0.7)"),
        opt("--admission_guard", dest="admission_guard", default='0.6,0.6,0.1,0.8,0.6', metavar="LIST", cmds=AP,
            help="shares of an added candidate's reads: tailed, in a TES cluster, internally primed, in a TSS "
                 "cluster outside genes, in a TSS cluster (default: 0.6,0.6,0.1,0.8,0.6)"),
        opt("--admission_min_datasets", dest="admission_min_datasets", type=int, default=2, metavar="N", cmds=AP,
            help="datasets in which an added candidate needs reads, in a run of several (default: 2)"),
    ]),
    ("Counting", [
        opt("--quant_requires_end_support", dest="quant_requires_end_support", action=Switch, default=True, cmds=AP,
            help="leave reads without a TSS or TES cluster out of the counts (default: on)"),
    ]),
    ("Small RNAs", [
        opt("--small_rna", dest="small_rna", action=Switch, default=None, cmds=A, main=A,
            help="count annotated small RNAs; needs the genome FASTA (default: on)"),
        opt("--small_rna_types", dest="small_rna_types",
            default='rRNA,snRNA,snoRNA,scaRNA,misc_RNA,vault_RNA,Y_RNA,Mt_tRNA,ribozyme,sRNA,scRNA',
            metavar="LIST", cmds=A,
            help="gene types of small RNAs, comma-separated (default: "
                 "rRNA,snRNA,snoRNA,scaRNA,misc_RNA,vault_RNA,Y_RNA,Mt_tRNA,ribozyme,sRNA,scRNA)"),
        opt("--small_rna_min_read_length", dest="small_rna_min_read_length", type=int, default=50, metavar="N",
            cmds=A,
            help="shortest read kept on a small RNA, in bases (default: 50)"),
        opt("--small_rna_max_length", dest="small_rna_max_length", type=int, default=500, metavar="N", cmds=A,
            help="longest transcript treated as a small RNA, in bases (default: 500)"),
    ]),
    ("Fusions", [
        opt("--fusions", dest="fusions", action=Switch, default=True, cmds=A, main=A,
            help="call fusion transcripts (default: on)"),
        opt("--fusion_spliced_joins", dest="fusion_spliced_joins", action=Switch, default=False, cmds=A,
            help="also cut a spliced alignment at an intron that joins two genes (experimental; default: off)"),
        opt("--fusion_rescue", dest="fusion_rescue", action=Switch, default=True, cmds=AP,
            help="place read parts the aligner left unaligned; needs sorted, indexed BAMs and the genome FASTA "
                 "(default: on)"),
        opt("--fusion_rescue_config", dest="fusion_rescue_config", metavar="CSV", cmds=AP,
            help="sample sheet of the run's BAMs, when they have moved"),
        opt("--fusion_min_reads", dest="fusion_min_reads", type=int, default=2, metavar="N", cmds=AP,
            help="reads a fusion needs with a breakpoint on an annotated exon boundary (default: 2)"),
        opt("--fusion_min_reads_intraexonic", dest="fusion_min_reads_intraexonic", type=int, default=4,
            metavar="N", cmds=AP,
            help="reads a fusion needs with no breakpoint on an exon boundary (default: 4)"),
        opt("--fusion_min_reads_unannotated", dest="fusion_min_reads_unannotated", type=int, default=5,
            metavar="N", cmds=AP,
            help="reads a fusion with an unannotated partner needs (default: 5)"),
        opt("--fusion_single_dataset_min_reads", dest="fusion_single_dataset_min_reads", type=int, default=5,
            metavar="N", cmds=AP,
            help="reads a fusion found in one dataset of several needs (default: 5)"),
        opt("--fusion_min_mapq", dest="fusion_min_mapq", type=int, default=20, metavar="N", cmds=AP,
            help="mapping quality each part of a fusion read needs (default: 20)"),
        opt("--fusion_min_fraction", dest="fusion_min_fraction", type=float, default=0.01, metavar="X", cmds=AP,
            help="fusion reads as a fraction of the less expressed partner's reads (default: 0.01)"),
        opt("--fusion_max_partners", dest="fusion_max_partners", type=int, default=3, metavar="N", cmds=AP,
            help="a gene joined to more partners than this is listed as promiscuous (default: 3)"),
        opt("--fusion_hub_dominance", dest="fusion_hub_dominance", type=float, default=20.0, metavar="X", cmds=AP,
            help="a call with this many times the median reads of such a gene's other calls passes (default: 20)"),
        opt("--fusion_hub_dominance_min_reads", dest="fusion_hub_dominance_min_reads", type=int, default=10,
            metavar="N", cmds=AP,
            help="reads such a call needs (default: 10)"),
        opt("--fusion_foldback_min_reads", dest="fusion_foldback_min_reads", type=int, default=5, metavar="N",
            cmds=AP,
            help="reads a call between overlapping opposite-strand neighbours needs to be no fold-back (default: "
                 "5)"),
        opt("--fusion_foldback_gap", dest="fusion_foldback_gap", type=int, default=1000, metavar="N", cmds=AP,
            help="bp between the two parts beyond which such a call is no fold-back (default: 1000)"),
        opt("--fusion_readthrough_distance", dest="fusion_readthrough_distance", type=int, default=1000000,
            metavar="N", cmds=AP,
            help="a same-strand partner this close downstream is read-through, not a fusion (default: 1000000)"),
        opt("--fusion_unannotated_partners", dest="fusion_unannotated_partners", action=Switch, default=True,
            cmds=AP,
            help="list fusions with an antisense, intronic or intergenic partner (default: on)"),
        opt("--fusion_pseudogene_partners", dest="fusion_pseudogene_partners", action=Switch, default=True, cmds=AP,
            help="call fusions with a pseudogene partner (default: on)"),
    ]),
    ("Run", [
        opt("--threads", "-t", dest="threads", type=int, default=os.cpu_count() or 4, metavar="N", cmds=A, main=A,
            help="threads (default: all)"),
        opt("--batch_size", dest="batch_size", type=int, default=1000, metavar="N", cmds=A,
            help="reads per batch (default: 1000)"),
        opt("--verbosity", "-v", dest="verbosity", type=int, default=1, metavar="N", cmds=A, main=A,
            help="0 quiet, 1 normal, 2 debug (default: 1)"),
    ]),
]


# The sections of the short help; the options of the other groups are listed there under one heading.
SHORT_SECTIONS = ("Input", "Reference", "Output", "Alignment of reads", "Run")


class HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Help with the descriptions in one column, at most 100 characters wide, and one spelling per option."""

    def __init__(self, prog):
        super().__init__(prog, max_help_position=30, width=min(shutil.get_terminal_size().columns, 100))

    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        long = [o for o in action.option_strings if o.startswith("--")]
        if isinstance(action, Switch):
            return f"{long[0]}, --no-{long[0][2:]}"
        names = [o for o in action.option_strings if not o.startswith("--")] + long[:1]
        if action.nargs == 0:
            return ", ".join(names)
        value = self._format_args(action, self._get_default_metavar_for_optional(action))
        return ", ".join(f"{name} {value}" for name in names)


def parser_for(command, usage, description, epilog, show_all, version=None):
    """The parser of a command. Its help lists every option with `show_all`, else those of everyday use."""
    parser = argparse.ArgumentParser(prog=f"isoforge {command}", usage=usage, description=description, epilog=epilog,
                                     formatter_class=HelpFormatter, add_help=False)
    sections = {}
    for title, options in GROUPS:
        if not show_all and title not in SHORT_SECTIONS:
            title = "Settings"
        for flags, cmds, main, settings in options:
            if command not in cmds:
                continue
            if title not in sections:
                sections[title] = parser.add_argument_group(title)
            settings = dict(settings)
            if not show_all and command not in main:
                settings["help"] = argparse.SUPPRESS
            sections[title].add_argument(*flags, **settings)
    group = parser.add_argument_group("Help")
    group.add_argument("-h", "--help", action="help", help="show the main options")
    group.add_argument("--help-all", action="help", help="show every option")
    if version:
        group.add_argument("-V", "--version", action="version", version=version, help="show the version")
    return parser
