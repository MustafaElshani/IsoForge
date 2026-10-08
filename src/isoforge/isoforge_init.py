#!/usr/bin/env python3
"""isoforge index — build the annotation index from a GTF file.

    isoforge index --gtf annotation.gtf --genome hg38 --annot gencode_v49

The index is written to `isoforge_index/` beside the GTF unless `--o` names a directory. An index that exists for
the same genome build and annotation name is kept; `--force` builds it again. `isoforge annotate --annotation`
builds the index of a GTF by itself and uses it again in later runs.
"""

import argparse
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from isoforge.options import HelpFormatter

__version__ = "1.0.0"

_INDEX_DIR_NAME = "isoforge_index"


def get_args():
    parser = argparse.ArgumentParser(
        prog="isoforge index",
        usage="isoforge index --gtf GTF --genome NAME --annot NAME [--o DIR]",
        description="Build the annotation index that `isoforge annotate` reads. Needed once per annotation.",
        epilog="Example:\n  isoforge index --gtf gencode.v49.annotation.gtf --genome hg38 --annot gencode_v49\n",
        formatter_class=HelpFormatter,
        add_help=False,
    )
    required = parser.add_argument_group("Required")
    required.add_argument("--gtf", "--f", dest="gtf", required=True, metavar="GTF",
                          help="annotation GTF with gene, transcript and exon records")
    required.add_argument("--genome", "--g", dest="genome", required=True, metavar="NAME",
                          help="genome build name, e.g. hg38; given again to `isoforge annotate`")
    required.add_argument("--annot", "--a", dest="annot", required=True, metavar="NAME",
                          help="annotation name, e.g. gencode_v49")
    optional = parser.add_argument_group("Options")
    optional.add_argument("--o", dest="output_dir", default=None, metavar="DIR",
                          help=f"index directory (default: {_INDEX_DIR_NAME}/ beside the GTF)")
    optional.add_argument("--min_length", "--l", dest="min_length", type=int, default=0, metavar="N",
                          help="leave out transcripts shorter than this, in bp (default: 0)")
    optional.add_argument("--force", action="store_true", help="build again an index that exists")
    optional.add_argument("-h", "--help", action="help", help="show this help")
    return parser.parse_args()


def _write_metadata(out_dir: Path, genome_build: str, annot_name: str,
                    gtf_path: str, min_length: int, version: str):
    """Write metadata.json alongside the Parquet files for validation."""
    meta = {
        "isoforge_version": version,
        "genome_build": genome_build,
        "annot_name": annot_name,
        "gtf_source": str(gtf_path),
        "min_length": min_length,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    meta_path = out_dir / "metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return meta_path


def index_metadata(out_dir):
    """The metadata of the index in a directory; None when the directory holds no complete index."""
    out_dir = Path(out_dir)
    meta_path = out_dir / "metadata.json"
    if not meta_path.exists() or not (out_dir / "edges.parquet").exists():
        return None
    try:
        return json.loads(meta_path.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def _check_existing_index(out_dir: Path, genome_build: str, annot_name: str) -> bool:
    """Return True if an existing index matches genome_build + annot_name."""
    meta = index_metadata(out_dir)
    return bool(meta) and meta.get("genome_build") == genome_build and meta.get("annot_name") == annot_name


def build_index(gtf_path, genome_build, annot_name, out_dir, min_length=0):
    """Builds the index of a GTF in out_dir and writes its metadata."""
    try:
        import isoforge_core
    except ImportError:
        print(
            "ERROR: the isoforge_core Rust extension is not installed.\n"
            "  Install IsoForge from its repository: bash install.sh (or pip install .)",
            file=sys.stderr,
        )
        sys.exit(1)

    gtf_path, out_dir = Path(gtf_path), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"IsoForge v{isoforge_core.version()} — Building annotation index")
    print(f"  GTF:        {gtf_path}")
    print(f"  Genome:     {genome_build}")
    print(f"  Annotation: {annot_name}")
    print(f"  Min length: {min_length} bp")
    print(f"  Output:     {out_dir}")
    print()

    t0 = time.time()
    try:
        isoforge_core.init_from_gtf(
            gtf_path=str(gtf_path),
            genome_build=genome_build,
            annot_name=annot_name,
            output_dir=str(out_dir),
            min_length=min_length,
        )
    except Exception as e:
        print(f"ERROR: Index build failed: {e}", file=sys.stderr)
        sys.exit(1)

    _write_metadata(out_dir, genome_build, annot_name, str(gtf_path), min_length, isoforge_core.version())

    print(f"  Completed in {time.time() - t0:.1f}s")
    print(f"  Output files:")
    for f in sorted(out_dir.glob("*.parquet")):
        size = f.stat().st_size
        if size > 1_000_000:
            print(f"    {f.name:<40} {size/1_000_000:.1f} MB")
        else:
            print(f"    {f.name:<40} {size/1_000:.1f} KB")
    for f in sorted(out_dir.glob("*.json")):
        print(f"    {f.name:<40} {f.stat().st_size/1_000:.1f} KB")
    print()
    print(f"  ✓ Annotation index ready at: {out_dir}")


def index_for_gtf(gtf, genome_build, fallback_dir):
    """The index of a GTF, built when there is none yet.

    An index beside the GTF is used when it was built from this GTF and is not older than it: `isoforge_index/`,
    else `<GTF name>_isoforge_index/`. A new index is written to the second, or into fallback_dir when the GTF's
    directory cannot be written.
    """
    gtf_path = Path(gtf).resolve()
    if not gtf_path.exists():
        print(f"ERROR: annotation GTF not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)
    name = re.sub(r"\.(gtf|gff3?)(\.gz)?$", "", gtf_path.name, flags=re.IGNORECASE)

    def is_index_of_gtf(directory):
        meta = index_metadata(directory)
        return (meta is not None and Path(str(meta.get("gtf_source", ""))).name == gtf_path.name
                and (directory / "metadata.json").stat().st_mtime >= gtf_path.stat().st_mtime)

    own = gtf_path.parent / f"{name}_{_INDEX_DIR_NAME}"
    elsewhere = Path(fallback_dir).resolve() / f"{name}_{_INDEX_DIR_NAME}"
    for directory in (gtf_path.parent / _INDEX_DIR_NAME, own, elsewhere):
        if is_index_of_gtf(directory):
            print(f"[IsoForge] Index of {gtf_path.name}: {directory}")
            return directory
    out_dir = own if os.access(gtf_path.parent, os.W_OK) else elsewhere
    build_index(gtf_path, genome_build, name, out_dir)
    return out_dir


def main():
    args = get_args()

    # Validate GTF exists
    gtf_path = Path(args.gtf).resolve()
    if not gtf_path.exists():
        print(f"ERROR: GTF file not found: {gtf_path}", file=sys.stderr)
        sys.exit(1)

    # Default output: isoforge_index/ next to the GTF
    if args.output_dir:
        out_dir = Path(args.output_dir).resolve()
    else:
        out_dir = gtf_path.parent / _INDEX_DIR_NAME

    # Skip-if-exists (unless --force)
    if not args.force and _check_existing_index(out_dir, args.genome, args.annot):
        meta = json.loads((out_dir / "metadata.json").read_text())
        print(f"[IsoForge] Index already exists and matches:")
        print(f"   Directory:   {out_dir}")
        print(f"   Genome:      {meta['genome_build']}")
        print(f"   Annotation:  {meta['annot_name']}")
        print(f"   Created:     {meta.get('created_at', 'unknown')}")
        print(f"   Skipping rebuild. Use --force to rebuild.")
        return

    build_index(gtf_path, args.genome, args.annot, out_dir, args.min_length)
    print(f"  Next step: isoforge annotate --db {out_dir} ...")


if __name__ == "__main__":
    main()
