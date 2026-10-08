#!/usr/bin/env python3
"""isoforge export — write coverage tracks or the read table again from a finished run.

    isoforge export tracks     --read_annot X_isoforge_read_annot.parquet --fai genome.fa.fai --o X
    isoforge export read_annot --read_annot X_isoforge_read_annot.parquet --o X

Both are made from the read-annotation Parquet, which holds every read of the run. The tracks are the files
`isoforge annotate` writes. The read table holds the same reads and columns, with the A content at the precision
the Parquet stores.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from isoforge.options import HelpFormatter

__version__ = "1.0.0"

# The column order of the read annotation table, as written by `isoforge annotate`.
TSV_COLUMNS = [
    "read_name", "dataset", "genome_build", "chromosome", "read_start", "read_end", "strand",
    "n_exons", "read_length", "gene_ID", "transcript_ID", "annot_gene_id", "annot_transcript_id",
    "annot_gene_name", "annot_transcript_name", "gene_novelty", "transcript_novelty",
    "ISM_subtype", "fraction_As", "custom_label", "allelic_label", "start_support", "end_support",
    "tss_site", "tss_feature", "polya_site", "polya_feature", "polya_signal", "cluster_supported",
    "polya_length",
]


EPILOG = """\
Commands:
  tracks       read-end coverage tracks per dataset and strand (bedGraph; BigWig with bedGraphToBigWig)
  read_annot   the read table as TSV

Examples:
  isoforge export tracks --read_annot results/project_isoforge_read_annot.parquet \\
      --fai hg38.fa.fai --o results/project
  isoforge export read_annot --read_annot results/project_isoforge_read_annot.parquet \\
      --whitelist results/project_filter.csv --o results/project

`isoforge export <command> --help` shows the options of a command.
"""


def get_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="isoforge export",
        usage="isoforge export {tracks,read_annot} --read_annot PARQUET --o PREFIX [options]",
        description="Write coverage tracks or the read table again from the read-annotation Parquet of a finished run.",
        epilog=EPILOG,
        formatter_class=HelpFormatter,
        add_help=False,
    )
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)
    sub = parser.add_subparsers(dest="what", required=True, metavar="{tracks,read_annot}", help=argparse.SUPPRESS)

    t = sub.add_parser("tracks", prog="isoforge export tracks", formatter_class=HelpFormatter, add_help=False,
                       usage="isoforge export tracks --read_annot PARQUET --fai FAI --o PREFIX [--no-bigwig]",
                       description="Write the read-end coverage tracks per dataset and strand.")
    t = t.add_argument_group("Options")
    t.add_argument("--read_annot", required=True, metavar="PARQUET", help="<prefix>_isoforge_read_annot.parquet of the run")
    t.add_argument("--fai", required=True, metavar="FAI", help=".fai index of the genome FASTA, for the chromosome sizes")
    t.add_argument("--o", dest="outprefix", required=True, metavar="PREFIX",
                   help="output prefix; the files go to <prefix>_TSS_TES_Clusters/counts/")
    t.add_argument("--no-bigwig", dest="bigwig", action="store_false", default=True, help="write bedGraph only")
    t.add_argument("-h", "--help", action="help", help="show this help")

    r = sub.add_parser("read_annot", prog="isoforge export read_annot", formatter_class=HelpFormatter, add_help=False,
                       usage="isoforge export read_annot --read_annot PARQUET --o PREFIX [options]",
                       description="Write the read table as TSV.")
    r = r.add_argument_group("Options")
    r.add_argument("--read_annot", required=True, metavar="PARQUET", help="<prefix>_isoforge_read_annot.parquet of the run")
    r.add_argument("--o", dest="outprefix", required=True, metavar="PREFIX", help="output prefix")
    r.add_argument("--whitelist", default=None, metavar="CSV",
                   help="<prefix>_filter.csv: write only the reads of these transcripts")
    r.add_argument("--split_unsupported", action="store_true",
                   help="write reads that are not whole reads of their model to a second table")
    r.add_argument("--maxFracA", dest="max_frac_a", type=float, default=None, metavar="X",
                   help="with --split_unsupported: also move reads without a poly(A) signal whose A fraction "
                        "is above this")
    r.add_argument("-h", "--help", action="help", help="show this help")
    return parser.parse_args(argv)


def _load(path, columns=None):
    import pandas as pd
    # nullable types keep integer columns with missing values as integers
    return pd.read_parquet(path, columns=columns, dtype_backend="numpy_nullable")


def export_tracks(args):
    import pandas as pd

    df = _load(args.read_annot, columns=["dataset", "chromosome", "strand", "tss_site", "polya_site"])
    out_dir = Path(f"{args.outprefix}_TSS_TES_Clusters") / "counts"
    out_dir.mkdir(parents=True, exist_ok=True)

    sizes = [ln.split("\t")[:2] for ln in Path(args.fai).read_text().splitlines() if ln.strip()]
    chrom_sizes = out_dir / "chrom.sizes"
    chrom_sizes.write_text("".join(f"{c}\t{n}\n" for c, n in sizes))

    bgtobw = shutil.which("bedGraphToBigWig") if args.bigwig else None
    datasets = sorted(df.dataset.unique())
    written = 0
    for label, site in (("tss", "tss_site"), ("polyA", "polya_site")):
        sub_all = df[df[site].notna()]
        for ds in [*datasets, "all"]:
            sub = sub_all if ds == "all" else sub_all[sub_all.dataset == ds]
            for strand, tag in (("+", "pos"), ("-", "neg")):
                s = sub[sub.strand == strand]
                path = out_dir / f"{ds}_{label}_counts_{tag}.bedGraph"
                counts = (s.groupby(["chromosome", site]).size().reset_index(name="n")
                          .sort_values(["chromosome", site]))
                with open(path, "w") as fh:
                    for chrom, pos, n in counts.itertuples(index=False):
                        fh.write(f"{chrom}\t{int(pos) - 1}\t{int(pos)}\t{n}\n")
                written += 1
                if bgtobw and len(counts):
                    bw = str(path).replace(".bedGraph", ".bw")
                    if subprocess.run([bgtobw, str(path), str(chrom_sizes), bw],
                                      capture_output=True).returncode == 0:
                        written += 1
    print(f"[export] {written} files in {out_dir}"
          f"{'' if bgtobw else '  (bedGraphToBigWig not found, bedGraph only)'}")


def export_read_annot(args):
    import pandas as pd

    df = _load(args.read_annot)
    missing = [c for c in TSV_COLUMNS if c not in df.columns]
    if missing:
        sys.exit(f"ERROR: {args.read_annot} lacks columns: {', '.join(missing)}")
    df = df[TSV_COLUMNS]

    if args.whitelist:
        keep = set()
        for line in Path(args.whitelist).read_text().splitlines():
            parts = line.split(",")
            if len(parts) > 1 and parts[1].strip().lstrip("-").isdigit():
                keep.add(int(parts[1]))
        if keep:
            df = df[df.transcript_ID.isin(keep)]

    main = df
    unsupported = None
    if args.split_unsupported:
        # a read without a value is not a whole read
        failed = ~df.cluster_supported.fillna(False).astype(bool)
        if args.max_frac_a is not None:
            failed = failed | ((df.fraction_As > args.max_frac_a).fillna(False).astype(bool) & df.polya_signal.isna())
        main, unsupported = df[~failed], df[failed]

    def write(frame, path):
        frame.to_csv(path, sep="\t", index=False, na_rep="")
        print(f"[export] {len(frame):,} reads -> {path}")

    write(main, f"{args.outprefix}_isoforge_read_annot.tsv")
    if unsupported is not None:
        write(unsupported, f"{args.outprefix}_isoforge_read_annot_unsupported.tsv")


def main(argv=None):
    args = get_args(argv)
    if args.what == "tracks":
        export_tracks(args)
    else:
        export_read_annot(args)


if __name__ == "__main__":
    sys.exit(main())
