#!/usr/bin/env python3
"""Run variant calling locally for a task JSON/dict."""

import argparse
import json
import os
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from niome_subnet.genomics.variant_caller import run_variant_calling

try:
    from niome_subnet.genomics.compare import (
        DEFAULT_TRUTH_ANNOTATIONS,
        DEFAULT_TRUTH_VCF,
        compare_submission,
    )
except ImportError:
    compare_submission = None


def main() -> None:
    parser = argparse.ArgumentParser(description="Run CFTR variant calling for a task")
    parser.add_argument(
        "--task-json",
        help="Path to a JSON file containing the task payload",
    )
    parser.add_argument(
        "--output-vcf",
        default="output.vcf",
        help="Path to write the generated VCF",
    )
    parser.add_argument(
        "--output-annotations",
        default="output.annotations.json",
        help="Path to write CFTR annotations JSON",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Score output against data/truth.vcf and data/annotations.json",
    )
    parser.add_argument(
        "--ref",
        default=None,
        help="Reference FASTA used when --compare is set (default: data/ref.fa)",
    )
    args = parser.parse_args()

    if args.task_json:
        with open(args.task_json, encoding="utf-8") as handle:
            task_data = json.load(handle)
    else:
        task_data = json.load(sys.stdin)

    vcf_content, annotations = run_variant_calling(task_data)

    with open(args.output_vcf, "w", encoding="utf-8") as handle:
        handle.write(vcf_content)

    with open(args.output_annotations, "w", encoding="utf-8") as handle:
        json.dump(annotations, handle, indent=2)

    variant_lines = [line for line in vcf_content.splitlines() if line and not line.startswith("#")]
    print(f"Wrote {args.output_vcf} ({len(variant_lines)} variant records)")
    print(f"Wrote {args.output_annotations} ({len(annotations)} annotated variants)")

    if args.compare:
        if compare_submission is None:
            raise SystemExit("compare_submission is unavailable (missing dependencies)")
        result = compare_submission(
            miner_vcf_path=args.output_vcf,
            miner_annotations_path=args.output_annotations,
            truth_vcf_path=DEFAULT_TRUTH_VCF,
            truth_annotations_path=DEFAULT_TRUTH_ANNOTATIONS,
            ref_fasta_path=args.ref or os.path.join(PROJECT_ROOT, "data", "ref.fa"),
        )
        print(f"Compared against {DEFAULT_TRUTH_VCF} and {DEFAULT_TRUTH_ANNOTATIONS}")
        print(f"Final score: {result.final_score:.4f} (VCF {result.vcf_score:.4f}, annotations {result.annotation_score:.4f})")


if __name__ == "__main__":
    main()
