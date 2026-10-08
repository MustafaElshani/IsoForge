#!/usr/bin/env python3
"""Alignment of read files for `isoforge annotate`.

Reads given as FASTQ or FASTA are aligned to the genome with the minimap2 built into IsoForge, into one sorted,
indexed BAM per sample in `<prefix>_alignments/`.
"""

import os
import re
import time
from pathlib import Path

READ_SUFFIXES = re.compile(r"\.(fastq|fq|fasta|fa|fna)(\.gz|\.bgz)?$", re.IGNORECASE)


def is_read_file(path):
    """Whether a file name is that of a FASTQ or FASTA file."""
    return bool(READ_SUFFIXES.search(str(path)))


def sample_name(path):
    """The sample name of a file: its name without the directory and the suffix."""
    name = os.path.basename(str(path))
    if is_read_file(name):
        return READ_SUFFIXES.sub("", name)
    return Path(name).stem


def align_samples(args, samples, out_prefix):
    """Aligns the read file of each sample, given as (reads, name), and returns the BAM written for each."""
    import isoforge_core

    bam_dir = Path(f"{out_prefix}_alignments")
    bam_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(reads, str(bam_dir / (re.sub(r"[^\w.\-]+", "_", name) + ".bam"))) for reads, name in samples]
    threads = args.align_threads or args.threads
    print(f"[IsoForge] Aligning {len(jobs)} read file(s) to {args.genome_fasta} on {threads} threads", flush=True)
    t0 = time.time()
    # a saved minimap2 index is read; otherwise the index is built from the FASTA (and saved when asked)
    saved_index = args.mm_index and os.path.exists(args.mm_index)
    # with annotated junctions (--junc_bed) the splice-site scores are used by annotation only
    align_scores = None if args.junc_bed else args.splice_scores
    stats = isoforge_core.align_reads(
        args.mm_index if saved_index else args.genome_fasta, jobs,
        preset="splice" if args.platform == "ont" else "splice:hq", k=args.k,
        forward_only=args.oriented, secondary=args.secondary, md=True, splice_scores=align_scores,
        splice_score_scale=args.splice_score_scale, sort=not args.no_sort, tmp_dir=str(bam_dir),
        sort_mem_gb=args.sort_mem_gb, junc_bed=args.junc_bed, copy_comment=args.copy_comment,
        index_out=None if saved_index else args.mm_index, threads=threads)
    for (_, name), (n, unmapped, primary, supp, _) in zip(samples, stats):
        print(f"[IsoForge]   {name}: {n:,} reads, {primary:,} aligned, {supp:,} supplementary, {unmapped:,} unmapped")
    print(f"[IsoForge] Alignment done in {time.time() - t0:.1f}s", flush=True)
    return [bam for _, bam in jobs]


def remove_alignments(bams):
    """Removes BAMs written by align_samples with their indexes, and their directory when it is left empty."""
    for bam in bams:
        for path in (bam, bam + ".bai"):
            try:
                os.remove(path)
            except OSError:
                pass
    for directory in {os.path.dirname(bam) for bam in bams}:
        try:
            os.rmdir(directory)
        except OSError:
            pass
