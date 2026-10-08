#!/usr/bin/env python3
"""
CAGE peak utilities for IsoForge.

Downloads FANTOM5 CAGE peak BED files and validates CAGE peak files given by
the user (`isoforge annotate --cage_peaks`).

Usage (Python):
    from isoforge.cage_utils import download_fantom5_cage
    bed_path = download_fantom5_cage("hg38", "/path/to/output/")

Usage (CLI):
    python -m isoforge.cage_utils --genome hg38 --output-dir /path/to/output/

Supported genomes:
    - hg38 (GRCh38): FANTOM5 reprocessed CAGE peaks
"""

import gzip
import os
import shutil
import urllib.request
from pathlib import Path

# FANTOM5 CAGE peak URLs
FANTOM5_URLS = {
    "hg38": "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz",
    "hg19": "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg19_latest/extra/CAGE_peaks/hg19_fair+new_CAGE_peaks_phase1and2.bed.gz",
}


def download_fantom5_cage(
    genome_build: str = "hg38",
    output_dir: str = ".",
    keep_gz: bool = True,
) -> str:
    """
    Download FANTOM5 CAGE peaks for the specified genome build.

    Args:
        genome_build: 'hg38' or 'hg19'
        output_dir: Directory to save the file
        keep_gz: If True, keep the gzipped file (IsoForge reads .gz natively).
                 If False, decompress to plain BED.

    Returns:
        Path to the downloaded (or decompressed) file.
    """
    if genome_build not in FANTOM5_URLS:
        raise ValueError(
            f"Unsupported genome build: {genome_build}. "
            f"Supported: {list(FANTOM5_URLS.keys())}"
        )

    url = FANTOM5_URLS[genome_build]
    filename = url.split("/")[-1]
    output_path = Path(output_dir) / filename

    if output_path.exists():
        print(f"[CAGE] File already exists: {output_path}")
        return str(output_path)

    os.makedirs(output_dir, exist_ok=True)
    print(f"[CAGE] Downloading FANTOM5 CAGE peaks for {genome_build}...")
    print(f"[CAGE] URL: {url}")

    urllib.request.urlretrieve(url, str(output_path))
    size_mb = output_path.stat().st_size / 1e6
    print(f"[CAGE] Downloaded: {output_path} ({size_mb:.1f} MB)")

    if not keep_gz:
        plain_path = output_path.with_suffix("")  # remove .gz
        print(f"[CAGE] Decompressing to {plain_path}...")
        with gzip.open(output_path, "rb") as f_in:
            with open(plain_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        output_path.unlink()
        return str(plain_path)

    return str(output_path)


def validate_cage_bed(path: str) -> dict:
    """
    Validate a CAGE peak BED file and return summary statistics.

    Args:
        path: Path to BED or BED.gz file

    Returns:
        Dict with 'n_peaks', 'n_chroms', 'has_strand', 'sample_lines'
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"CAGE peak file not found: {path}")

    opener = gzip.open if str(path).endswith(".gz") else open

    n_peaks = 0
    chroms = set()
    has_strand = False
    sample_lines = []

    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            n_peaks += 1
            chroms.add(fields[0])
            if len(fields) > 5 and fields[5] in ("+", "-"):
                has_strand = True
            if len(sample_lines) < 3:
                sample_lines.append(line)

    return {
        "n_peaks": n_peaks,
        "n_chroms": len(chroms),
        "has_strand": has_strand,
        "sample_lines": sample_lines,
    }


def main():
    """Command-line entry point: python -m isoforge.cage_utils."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Download FANTOM5 CAGE peaks for IsoForge 5′ end validation"
    )
    parser.add_argument(
        "--genome", default="hg38", choices=["hg38", "hg19"],
        help="Genome build (default: hg38)"
    )
    parser.add_argument(
        "--output-dir", default=".",
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "--decompress", action="store_true",
        help="Decompress to plain BED (default: keep .gz)"
    )

    args = parser.parse_args()
    path = download_fantom5_cage(args.genome, args.output_dir, keep_gz=not args.decompress)
    
    # Validate
    info = validate_cage_bed(path)
    print(f"[CAGE] Validation: {info['n_peaks']:,} peaks across {info['n_chroms']} chromosomes")
    print(f"[CAGE] Stranded: {'yes' if info['has_strand'] else 'no'}")
    print(f"[CAGE] Ready for IsoForge: cage_peaks_file='{path}'")


if __name__ == "__main__":
    main()
