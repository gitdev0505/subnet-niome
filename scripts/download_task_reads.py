#!/usr/bin/env python3
"""Download task FASTQs to data/read_1.fq and data/read_2.fq for offline runs."""

import argparse
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from niome_subnet.genomics.variant_caller import ensure_read_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download paired FASTQs from presigned URLs into data/",
    )
    parser.add_argument("--read1-url", required=True, help="Presigned URL for read 1")
    parser.add_argument("--read2-url", required=True, help="Presigned URL for read 2")
    parser.add_argument(
        "--dest-dir",
        default=os.path.join(PROJECT_ROOT, "data"),
        help="Output directory (default: data/)",
    )
    args = parser.parse_args()

    os.makedirs(args.dest_dir, exist_ok=True)
    read1_path = os.path.join(args.dest_dir, "read_1.fq")
    read2_path = os.path.join(args.dest_dir, "read_2.fq")

    for path in (read1_path, read2_path):
        if os.path.exists(path):
            os.remove(path)

    ensure_read_file(args.read1_url, read1_path)
    ensure_read_file(args.read2_url, read2_path)

    print(f"Wrote {read1_path}")
    print(f"Wrote {read2_path}")
    print("Update sample_task.json to use data/read_1.fq and data/read_2.fq, then run run_variant_calling.py")


if __name__ == "__main__":
    main()
