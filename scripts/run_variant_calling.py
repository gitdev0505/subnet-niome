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


if __name__ == "__main__":
    main()
