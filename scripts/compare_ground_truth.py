#!/usr/bin/env python3
"""Compare a miner VCF against validator ground truth in data/."""

import argparse
import json
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from niome_subnet.genomics.compare import (
    DEFAULT_READ1,
    DEFAULT_READ2,
    DEFAULT_REF_FASTA,
    DEFAULT_TRUTH_ANNOTATIONS,
    DEFAULT_TRUTH_VCF,
    compare_submission,
    map_vcf_to_cftr2_annotations,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Score miner output against data/truth.vcf and data/annotations.json",
    )
    parser.add_argument(
        "--miner-vcf",
        required=True,
        help="Path to the miner-generated VCF",
    )
    parser.add_argument(
        "--miner-annotations",
        help="Path to miner CFTR2 annotations JSON (optional)",
    )
    parser.add_argument(
        "--truth-vcf",
        default=DEFAULT_TRUTH_VCF,
        help=f"Validator truth VCF (default: {DEFAULT_TRUTH_VCF})",
    )
    parser.add_argument(
        "--truth-annotations",
        default=DEFAULT_TRUTH_ANNOTATIONS,
        help=f"Validator annotations JSON (default: {DEFAULT_TRUTH_ANNOTATIONS})",
    )
    parser.add_argument(
        "--ref",
        default=DEFAULT_REF_FASTA,
        help=f"Reference FASTA for normalization (default: {DEFAULT_REF_FASTA})",
    )
    parser.add_argument(
        "--read1",
        default=DEFAULT_READ1,
        help=f"Read 1 FASTQ for depth weighting (default: {DEFAULT_READ1})",
    )
    parser.add_argument(
        "--read2",
        default=DEFAULT_READ2,
        help=f"Read 2 FASTQ for depth weighting (default: {DEFAULT_READ2})",
    )
    parser.add_argument(
        "--response-time",
        type=float,
        help="Optional miner response time in seconds",
    )
    parser.add_argument(
        "--write-annotations",
        help="Write CFTR2-keyed annotations derived from the miner VCF to this path",
    )
    args = parser.parse_args()

    for label, path in (
        ("truth VCF", args.truth_vcf),
        ("truth annotations", args.truth_annotations),
        ("reference FASTA", args.ref),
        ("miner VCF", args.miner_vcf),
    ):
        if not os.path.exists(path):
            raise SystemExit(f"Missing {label}: {path}")

    for label, path in (("read 1", args.read1), ("read 2", args.read2)):
        if not os.path.exists(path):
            print(f"Warning: missing {label} FASTQ: {path} (depth weighting disabled)")

    if args.miner_annotations and not os.path.exists(args.miner_annotations):
        raise SystemExit(f"Missing miner annotations: {args.miner_annotations}")

    if args.write_annotations:
        with open(args.miner_vcf, encoding="utf-8") as handle:
            vcf_content = handle.read()
        derived = map_vcf_to_cftr2_annotations(
            vcf_content,
            truth_vcf_path=args.truth_vcf,
            truth_annotations_path=args.truth_annotations,
        )
        with open(args.write_annotations, "w", encoding="utf-8") as handle:
            json.dump(derived, handle, indent=2)
        print(f"Wrote {args.write_annotations} ({len(derived)} variants)")

    result = compare_submission(
        miner_vcf_path=args.miner_vcf,
        miner_annotations_path=args.miner_annotations,
        truth_vcf_path=args.truth_vcf,
        truth_annotations_path=args.truth_annotations,
        ref_fasta_path=args.ref,
        read1_path=args.read1,
        read2_path=args.read2,
        response_time=args.response_time,
        build_annotations=args.miner_annotations is None,
    )

    print(f"UID:               {result.uid}")
    print(f"Precision:         {result.precision:.4f}")
    print(f"Recall:            {result.recall:.4f}")
    print(f"F1:                {result.f1_score:.4f}")
    print(f"VCF score:         {result.vcf_score:.4f}")
    print(f"Annotation score:  {result.annotation_score:.4f}")
    print(f"Final score:       {result.final_score:.4f}")
    if result.response_time is not None:
        print(f"Response time:     {result.response_time:.2f}s")


if __name__ == "__main__":
    main()
