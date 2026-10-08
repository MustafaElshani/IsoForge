#!/usr/bin/env python3
"""The `isoforge` command.

    isoforge annotate     from FASTQ or sorted BAM: transcript models, ends, counts, fusions
    isoforge index        build the annotation index ahead of a run
    isoforge post         filter and write the outputs again from a kept run database
    isoforge export       write coverage tracks or the read table again from a finished run
"""

import sys

from isoforge.isoforge_main import __version__


def main():
    # Dispatch to the subcommand; each has its own argument parser.
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        _print_help()
        sys.exit(0 if len(sys.argv) >= 2 else 2)

    if sys.argv[1] in ("-V", "--version"):
        try:
            import isoforge_core
            print(f"IsoForge {isoforge_core.version()}")
        except ImportError:
            print(f"IsoForge {__version__} (isoforge_core not installed)")
        sys.exit(0)

    subcmd = sys.argv[1]

    # Remove the subcommand from argv: "isoforge annotate --db ..." becomes
    # prog="isoforge annotate", argv=["--db", ...].
    sys.argv = [f"isoforge {subcmd}"] + sys.argv[2:]

    if subcmd == "index":
        from isoforge.isoforge_init import main as _init_main
        _init_main()
    elif subcmd == "annotate":
        from isoforge.isoforge_main import main as _annotate_main
        _annotate_main()
    elif subcmd == "post":
        from isoforge.isoforge_post import main as _post_main
        _post_main()
    elif subcmd == "export":
        from isoforge.isoforge_export import main as _export_main
        _export_main()
    else:
        print(f"isoforge: unknown command '{subcmd}'", file=sys.stderr)
        print(file=sys.stderr)
        _print_help(file=sys.stderr)
        sys.exit(2)


def _print_help(file=None):
    file = file or sys.stdout
    try:
        import isoforge_core
        ver = isoforge_core.version()
    except ImportError:
        ver = __version__

    print(f"""\
IsoForge v{ver} — 5′ Isoforms · Ends · Counts · Fusions 3′

usage: isoforge <command> [options]

Commands:
  annotate     from FASTQ or sorted BAM: transcript models with their ends, counts, fusions
  index        build the annotation index ahead of a run (annotate builds it when needed)
  post         filter and write the outputs again from a kept run database
  export       write coverage tracks or the read table again from a finished run

  -h, --help      show this help
  -V, --version   show the version

Example:
  isoforge annotate --reads sample.fastq.gz --genome hg38.fa --annotation gencode.v49.gtf --o results/sample

`isoforge <command> --help` shows the options of a command.\
""", file=file)


if __name__ == "__main__":
    main()
