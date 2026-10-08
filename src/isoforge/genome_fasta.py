#!/usr/bin/env python3
"""The index of the genome FASTA.

IsoForge reads the genome through its `.fai` index. `ensure_fai` writes the index beside a FASTA that has none, in
the format of `samtools faidx`: name, length, offset of the first base, bases per line, bytes per line.
"""

import os
import sys


def ensure_fai(fasta):
    """Writes `<fasta>.fai` when it is missing. Stops with a message when the FASTA cannot be indexed."""
    fai = fasta + ".fai"
    if os.path.exists(fai):
        return
    if fasta.lower().endswith(".gz"):
        sys.exit(f"ERROR: {fasta} is compressed and has no index; index it with `samtools faidx` (bgzip only), "
                 "or give the FASTA uncompressed.")
    print(f"[IsoForge] Indexing the genome FASTA: {fai}", flush=True)
    records = []
    name = None
    offset = 0
    with open(fasta, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    records.append((name, length, start, line_bases, line_bytes))
                name = line[1:].split()[0].decode() if line[1:].split() else ""
                length, start, line_bases, line_bytes, short_line = 0, offset + len(line), 0, 0, False
            elif name is not None:
                bases = len(line.rstrip(b"\r\n"))
                if not line_bases and bases:
                    line_bases, line_bytes = bases, len(line)
                elif bases != line_bases or len(line) != line_bytes:
                    # only the last line of a sequence may differ from the others
                    if bases > line_bases or (short_line and bases):
                        sys.exit(f"ERROR: {fasta}: the lines of sequence {name} differ in length; "
                                 "the FASTA cannot be indexed.")
                    short_line = True
                elif short_line:
                    sys.exit(f"ERROR: {fasta}: the lines of sequence {name} differ in length; "
                             "the FASTA cannot be indexed.")
                length += bases
            offset += len(line)
    if name is not None:
        records.append((name, length, start, line_bases, line_bytes))
    if not records:
        sys.exit(f"ERROR: {fasta} holds no sequence.")
    try:
        with open(fai, "w") as out:
            for record in records:
                out.write("\t".join(str(x) for x in record) + "\n")
    except OSError as error:
        sys.exit(f"ERROR: cannot write {fai} ({error.strerror}); index the FASTA with `samtools faidx`, "
                 "or copy it to a directory you can write to.")
