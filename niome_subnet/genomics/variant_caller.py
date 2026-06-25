"""CFTR variant calling pipeline for miner tasks."""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import urllib.request
from typing import Any, Dict, Optional, Tuple

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_REF_PATH = os.path.join(PROJECT_ROOT, "data", "ref.fa")
UCSC_SEQUENCE_API = "https://api.genome.ucsc.edu/getData/sequence"
SAMPLE_NAME = "SAMPLE"
_CFTR_DRUGS = (
    "ivacaftor",
    "tezacaftor_ivacaftor",
    "elexacaftor_tezacaftor_ivacaftor",
    "lumacaftor_ivacaftor",
)


class VariantCallingError(RuntimeError):
    """Raised when variant calling cannot be completed."""


def parse_region(region: str) -> Tuple[str, int, int]:
    """Parse a region string like ``chr7:117480000-117670000``."""
    match = re.match(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$", region)
    if not match:
        raise VariantCallingError(f"Invalid genome region: {region}")
    return match.group("chrom"), int(match.group("start")), int(match.group("end"))


def _require_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise VariantCallingError(
            f"Required tool '{name}' not found in PATH. "
            "Install bwa, samtools, and bcftools (see docs/miner_guide.md)."
        )
    return path


def _run(cmd: list[str], *, shell: bool = False) -> None:
    result = subprocess.run(
        cmd,
        shell=shell,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if result.returncode != 0:
        raise VariantCallingError(
            f"Command failed ({result.returncode}): {' '.join(cmd)}\n{result.stderr.strip()}"
        )


def download_file(url: str, dest: str) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    urllib.request.urlretrieve(url, dest)


def fetch_reference_region(chrom: str, start: int, end: int, dest: str) -> None:
    """Download GRCh38/hg38 sequence for a region from UCSC."""
    url = f"{UCSC_SEQUENCE_API}?genome=hg38;chrom={chrom};start={start};end={end}"
    with urllib.request.urlopen(url, timeout=120) as response:
        payload = json.loads(response.read().decode())
    sequence = payload.get("dna")
    if not sequence:
        raise VariantCallingError(f"UCSC returned no sequence for {chrom}:{start}-{end}")

    os.makedirs(os.path.dirname(dest), exist_ok=True)
    with open(dest, "w", encoding="utf-8") as handle:
        handle.write(f">{chrom}\n")
        for offset in range(0, len(sequence), 80):
            handle.write(sequence[offset : offset + 80] + "\n")


def ensure_reference(
    chrom: str,
    start: int,
    end: int,
    ref_fasta: Optional[str] = None,
) -> Tuple[str, int]:
    """Return reference FASTA path and coordinate offset for regional fallbacks."""
    ref_path = ref_fasta or os.environ.get("NIOME_REF_FASTA", DEFAULT_REF_PATH)
    if os.path.exists(ref_path):
        return ref_path, 0

    fetch_reference_region(chrom, start, end, ref_path)
    return ref_path, start - 1


def index_reference(ref_fasta: str) -> None:
    if os.path.exists(ref_fasta + ".bwt"):
        return
    _run([_require_tool("bwa"), "index", ref_fasta])


def align_reads(ref_fasta: str, read1: str, read2: str, bam_path: str) -> None:
    bwa = _require_tool("bwa")
    samtools = _require_tool("samtools")
    os.makedirs(os.path.dirname(bam_path), exist_ok=True)
    read_group = f"@RG\\tID:{SAMPLE_NAME}\\tSM:{SAMPLE_NAME}\\tPL:ILLUMINA"
    align_cmd = (
        f'"{bwa}" mem -t 4 -R "{read_group}" "{ref_fasta}" "{read1}" "{read2}" '
        f'| "{samtools}" sort -o "{bam_path}"'
    )
    _run(align_cmd, shell=True)
    _run([samtools, "index", bam_path])


def call_variants(
    ref_fasta: str,
    bam_path: str,
    region: str,
    vcf_path: str,
) -> None:
    bcftools = _require_tool("bcftools")
    mpileup_cmd = (
        f'"{bcftools}" mpileup -f "{ref_fasta}" -r {region} "{bam_path}" '
        f'| "{bcftools}" call -mv -Ov -o "{vcf_path}"'
    )
    _run(mpileup_cmd, shell=True)


def vcf_to_text(vcf_path: str) -> str:
    with open(vcf_path, encoding="utf-8") as handle:
        return handle.read()


def adjust_vcf_coordinates(vcf_content: str, offset: int) -> str:
    """Shift VCF positions from a regional reference to absolute genomic coordinates."""
    if offset == 0:
        return vcf_content

    adjusted: list[str] = []
    for line in vcf_content.splitlines():
        if not line or line.startswith("#"):
            adjusted.append(line)
            continue
        fields = line.split("\t")
        fields[1] = str(int(fields[1]) + offset)
        adjusted.append("\t".join(fields))
    return "\n".join(adjusted) + ("\n" if adjusted and adjusted[-1] else "")


def _vcf_column_names(header_line: str) -> set[str]:
    """Return normalized VCF column names from a ``#CHROM`` header line."""
    if not header_line or not header_line.startswith("#CHROM"):
        return set()
    columns = header_line.split("\t")
    columns[0] = columns[0].lstrip("#")
    return set(columns)


def normalize_vcf_sample_column(
    vcf_content: str, sample_name: str = SAMPLE_NAME
) -> str:
    """Rename the genotype sample column to ``sample_name`` if needed."""
    lines = vcf_content.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("#CHROM"):
            continue
        fields = line.split("\t")
        if len(fields) >= 10 and fields[-1] != sample_name:
            fields[-1] = sample_name
            lines[index] = "\t".join(fields)
        break
    else:
        return vcf_content

    body = "\n".join(lines)
    return body + ("\n" if body and not body.endswith("\n") else "")


def build_empty_vcf(region: str, sample_name: str = SAMPLE_NAME) -> str:
    chrom, start, end = parse_region(region)
    return (
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID={chrom},length={end}>\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
    )


def _default_drug_response() -> Dict[str, str]:
    return {drug: "unknown" for drug in _CFTR_DRUGS}


def annotate_cftr_variants(vcf_content: str) -> Dict[str, Any]:
    """Build CFTR2-style annotations for variants present in the VCF."""
    annotations: Dict[str, Any] = {}
    for line in vcf_content.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue

        chrom, pos, variant_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        rsid = variant_id if variant_id.startswith("rs") else f"{chrom}:{pos}"
        annotations[rsid] = {
            "hgvs": f"{chrom}:g.{pos}{ref}>{alt.split(',')[0]}",
            "clinical_significance": "unknown",
            "drug_response": _default_drug_response(),
        }
    return annotations


def run_variant_calling(
    task_data: dict,
    work_dir: Optional[str] = None,
    ref_fasta: Optional[str] = None,
) -> Tuple[str, Dict[str, Any]]:
    """
    Run the CFTR variant-calling pipeline for a miner task.

    Returns:
        Tuple of (vcf_content, cftr_annotations).
    """
    task_id = task_data["task_id"]
    region = task_data["genome_context"]["region"]
    chrom, start, end = parse_region(region)

    work_dir = work_dir or os.path.join(PROJECT_ROOT, "data", "tasks", task_id)
    os.makedirs(work_dir, exist_ok=True)

    read1_path = os.path.join(work_dir, "read_1.fq")
    read2_path = os.path.join(work_dir, "read_2.fq")
    bam_path = os.path.join(work_dir, "aligned.bam")
    vcf_path = os.path.join(work_dir, "variants.vcf")

    task_input = task_data["input"]
    if not os.path.exists(read1_path):
        download_file(task_input["read1_fastq"], read1_path)
    if not os.path.exists(read2_path):
        download_file(task_input["read2_fastq"], read2_path)

    ref_path, coordinate_offset = ensure_reference(chrom, start, end, ref_fasta=ref_fasta)
    call_region = region if coordinate_offset == 0 else f"{chrom}:1-{end - start}"

    index_reference(ref_path)
    align_reads(ref_path, read1_path, read2_path, bam_path)
    call_variants(ref_path, bam_path, call_region, vcf_path)

    try:
        vcf_content = adjust_vcf_coordinates(vcf_to_text(vcf_path), coordinate_offset)
    except (subprocess.CalledProcessError, VariantCallingError):
        vcf_content = build_empty_vcf(region)

    vcf_content = normalize_vcf_sample_column(vcf_content, SAMPLE_NAME)

    # Ensure required VCF columns are present in the header.
    required_fields = set(task_data.get("output_spec", {}).get("required_fields", []))
    header_line = next(
        (line for line in vcf_content.splitlines() if line.startswith("#CHROM")), ""
    )
    if not header_line:
        raise VariantCallingError("VCF missing #CHROM header line")

    missing = required_fields - _vcf_column_names(header_line)
    if missing:
        raise VariantCallingError(f"VCF missing required columns: {sorted(missing)}")

    cftr_annotations = annotate_cftr_variants(vcf_content)
    return vcf_content, cftr_annotations
