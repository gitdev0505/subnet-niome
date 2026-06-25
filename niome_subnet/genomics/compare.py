"""Compare miner submissions against local validator ground truth files."""

from __future__ import annotations

import json
import os
import re
from typing import Any, Dict, Optional, Tuple

from niome_subnet.genomics.model import GroundTruth, MinerScore, MinerSubmission
from niome_subnet.genomics.scoring import create_mapping_file, score

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_DATA_DIR = os.path.join(PROJECT_ROOT, "data")
DEFAULT_TRUTH_VCF = os.path.join(DEFAULT_DATA_DIR, "truth.vcf")
DEFAULT_TRUTH_ANNOTATIONS = os.path.join(DEFAULT_DATA_DIR, "annotations.json")
DEFAULT_REF_FASTA = os.path.join(DEFAULT_DATA_DIR, "ref.fa")
DEFAULT_READ1 = os.path.join(DEFAULT_DATA_DIR, "read_1.fq")
DEFAULT_READ2 = os.path.join(DEFAULT_DATA_DIR, "read_2.fq")
CHR7_ACCESSION = "NC_000007.14"
_POSITION_TOLERANCE = 5


def _default_drug_response() -> Dict[str, str]:
    return {
        "ivacaftor": "unknown",
        "tezacaftor_ivacaftor": "unknown",
        "elexacaftor_tezacaftor_ivacaftor": "unknown",
        "lumacaftor_ivacaftor": "unknown",
    }


def _hgvs_matches_variant(hgvs: str, pos: int, ref: str, alt: str) -> bool:
    """Return whether a CFTR2 HGVS string describes a VCF variant."""
    suffix = hgvs.split(":")[-1] if ":" in hgvs else hgvs
    if not suffix.startswith("g."):
        return False

    body = suffix[2:]
    alt = alt.split(",")[0]

    snv = re.match(r"^(\d+)([ACGT]+)>([ACGT]+)$", body)
    if snv:
        return int(snv.group(1)) == pos and snv.group(2) == ref and snv.group(3) == alt

    range_del = re.match(r"^(\d+)_(\d+)del$", body)
    if range_del:
        start, end = int(range_del.group(1)), int(range_del.group(2))
        return start - _POSITION_TOLERANCE <= pos <= end + _POSITION_TOLERANCE

    single_del = re.match(r"^(\d+)del$", body)
    if single_del:
        del_pos = int(single_del.group(1))
        return abs(pos - del_pos) <= _POSITION_TOLERANCE

    return False


def _match_cftr2_id(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    truth_annotations: Dict[str, Any],
) -> Optional[str]:
    for cftr2_id, entry in truth_annotations.items():
        if _hgvs_matches_variant(entry.get("hgvs", ""), pos, ref, alt):
            return cftr2_id
    return None


def build_truth_variant_index(
    truth_vcf_path: str,
    truth_annotations_path: str,
) -> Dict[Tuple[str, int, str, str], str]:
    """Map truth VCF records to CFTR2 annotation IDs."""
    if not os.path.exists(truth_vcf_path) or not os.path.exists(truth_annotations_path):
        return {}

    with open(truth_annotations_path, encoding="utf-8") as handle:
        truth_annotations = json.load(handle)

    lookup: Dict[Tuple[str, int, str, str], str] = {}
    with open(truth_vcf_path, encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            cftr2_id = _match_cftr2_id(chrom, pos, ref, alt, truth_annotations)
            if cftr2_id:
                lookup[(chrom, pos, ref, alt)] = cftr2_id
    return lookup


def _vcf_to_hgvs(chrom: str, pos: str, ref: str, alt: str) -> str:
    accession = CHR7_ACCESSION if chrom == "chr7" else chrom
    return f"{accession}:g.{pos}{ref}>{alt.split(',')[0]}"


def map_vcf_to_cftr2_annotations(
    vcf_content: str,
    truth_vcf_path: str = DEFAULT_TRUTH_VCF,
    truth_annotations_path: str = DEFAULT_TRUTH_ANNOTATIONS,
) -> Dict[str, Any]:
    """Build miner-style annotations keyed by CFTR2 IDs from a VCF."""
    if not os.path.exists(truth_annotations_path):
        return {}

    with open(truth_annotations_path, encoding="utf-8") as handle:
        truth_annotations = json.load(handle)

    truth_index = build_truth_variant_index(truth_vcf_path, truth_annotations_path)
    annotations: Dict[str, Any] = {}

    for line in vcf_content.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue

        chrom, pos, variant_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
        pos_int = int(pos)
        cftr2_id = truth_index.get((chrom, pos_int, ref, alt))
        if cftr2_id is None:
            cftr2_id = _match_cftr2_id(chrom, pos_int, ref, alt, truth_annotations)
        if cftr2_id is None:
            continue

        truth_entry = truth_annotations.get(cftr2_id, {})
        annotations[cftr2_id] = {
            "hgvs": truth_entry.get("hgvs") or _vcf_to_hgvs(chrom, pos, ref, alt),
            "clinical_significance": "unknown",
            "drug_response": _default_drug_response(),
        }

    return annotations


def compare_submission(
    miner_vcf_path: str,
    miner_annotations_path: Optional[str] = None,
    truth_vcf_path: str = DEFAULT_TRUTH_VCF,
    truth_annotations_path: str = DEFAULT_TRUTH_ANNOTATIONS,
    ref_fasta_path: str = DEFAULT_REF_FASTA,
    read1_path: Optional[str] = DEFAULT_READ1,
    read2_path: Optional[str] = DEFAULT_READ2,
    response_time: Optional[float] = None,
    uid: int = 0,
    build_annotations: bool = True,
) -> MinerScore:
    """
    Score a miner VCF (and optional annotations) against local ground truth.

    Ground truth defaults to ``data/truth.vcf`` and ``data/annotations.json``.
    """
    with open(miner_vcf_path, encoding="utf-8") as handle:
        vcf_content = handle.read()

    miner_annotations: Optional[Dict[str, Any]] = None
    if miner_annotations_path and os.path.exists(miner_annotations_path):
        with open(miner_annotations_path, encoding="utf-8") as handle:
            miner_annotations = json.load(handle)
    elif build_annotations:
        miner_annotations = map_vcf_to_cftr2_annotations(
            vcf_content,
            truth_vcf_path=truth_vcf_path,
            truth_annotations_path=truth_annotations_path,
        )

    ground_truth = GroundTruth(
        truth_vcf=truth_vcf_path,
        ref=ref_fasta_path,
        cftr2_annotations=truth_annotations_path,
    )

    bam = ""
    if (
        read1_path
        and read2_path
        and os.path.exists(read1_path)
        and os.path.exists(read2_path)
        and os.path.exists(ref_fasta_path)
    ):
        bam = create_mapping_file(ref_fasta_path, read1_path, read2_path)

    return score(
        MinerSubmission(
            uid=uid,
            vcf_content=vcf_content,
            response_time=response_time,
            cftr_annotations=miner_annotations,
        ),
        ground_truth,
        bam,
    )
