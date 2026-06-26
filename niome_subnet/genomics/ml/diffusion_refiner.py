"""Diffusion-inspired genotype refiner for low-coverage CFTR variant calling."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from niome_subnet.genomics.ml.pileup import PileupSite, collect_pileup_sites

DEFAULT_STEPS = 20
DEFAULT_MIN_CONFIDENCE = 0.55
DEFAULT_MIN_ALT_FRACTION = 0.12
SAMPLE_NAME = "SAMPLE"


@dataclass
class RefinedCall:
    chrom: str
    pos: int
    ref: str
    alt: str
    genotype: str
    confidence: float
    depth: int


def _gt_probs_from_string(gt: str) -> np.ndarray:
    """Map VCF GT to [P(0/0), P(0/1), P(1/1)]."""
    mapping = {
        "0/0": np.array([0.95, 0.04, 0.01]),
        "0|0": np.array([0.95, 0.04, 0.01]),
        "0/1": np.array([0.05, 0.90, 0.05]),
        "0|1": np.array([0.05, 0.90, 0.05]),
        "1/0": np.array([0.05, 0.90, 0.05]),
        "1|0": np.array([0.05, 0.90, 0.05]),
        "1/1": np.array([0.01, 0.04, 0.95]),
        "1|1": np.array([0.01, 0.04, 0.95]),
        "./.": np.array([0.34, 0.33, 0.33]),
        ".|.": np.array([0.34, 0.33, 0.33]),
    }
    return mapping.get(gt, np.array([0.34, 0.33, 0.33]))


def _gt_string_from_probs(probs: np.ndarray) -> Tuple[str, float]:
    labels = ("0/0", "0/1", "1/1")
    index = int(np.argmax(probs))
    return labels[index], float(probs[index])


class DiffusionGenotypeRefiner:
    """Iteratively denoise genotype posteriors conditioned on pileup alt fraction."""

    def __init__(
        self,
        steps: int = DEFAULT_STEPS,
        min_confidence: float = DEFAULT_MIN_CONFIDENCE,
        min_alt_fraction: float = DEFAULT_MIN_ALT_FRACTION,
    ) -> None:
        self.steps = steps
        self.min_confidence = min_confidence
        self.min_alt_fraction = min_alt_fraction

    def denoise(self, gt_probs: np.ndarray, alt_fraction: float, depth: int) -> np.ndarray:
        """Reverse diffusion steps pulling genotype mass toward pileup likelihood."""
        x = np.asarray(gt_probs, dtype=float)
        if x.sum() <= 0:
            x = np.array([1 / 3, 1 / 3, 1 / 3])
        else:
            x = x / x.sum()

        depth_weight = min(depth / 14.0, 1.0)
        for step in range(self.steps, 0, -1):
            temperature = step / self.steps
            logits = np.log(x + 1e-9)
            for genotype_index, expected_af in enumerate((0.0, 0.5, 1.0)):
                residual = alt_fraction - expected_af
                logits[genotype_index] -= (
                    (1.0 - temperature) * depth_weight * (residual ** 2) * depth / 10.0
                )
            x = np.exp(logits - logits.max())
            x /= x.sum()
        return x

    def refine_site(
        self,
        site: PileupSite,
        prior_gt: Optional[str] = None,
    ) -> Optional[RefinedCall]:
        if site.alt_allele is None or site.alt_fraction < self.min_alt_fraction:
            return None

        prior = _gt_probs_from_string(prior_gt) if prior_gt else np.array([0.7, 0.25, 0.05])
        posterior = self.denoise(prior, site.alt_fraction, site.depth)
        genotype, confidence = _gt_string_from_probs(posterior)

        if genotype == "0/0" or confidence < self.min_confidence:
            return None

        return RefinedCall(
            chrom=site.chrom,
            pos=site.pos,
            ref=site.ref,
            alt=site.alt_allele,
            genotype=genotype,
            confidence=confidence,
            depth=site.depth,
        )

    def refine_calls(
        self,
        pileups: List[PileupSite],
        existing: Dict[Tuple[int, str, str], str],
    ) -> List[RefinedCall]:
        """Refine bcftools calls and recover missed variants from pileups."""
        refined: List[RefinedCall] = []

        for site in pileups:
            if site.alt_allele is None:
                continue
            key = (site.pos, site.ref, site.alt_allele)
            call = self.refine_site(site, existing.get(key))
            if call is None:
                continue
            refined.append(call)

        return refined


def _parse_vcf_records(
    vcf_content: str,
) -> Tuple[List[str], Dict[Tuple[int, str, str], str]]:
    header: List[str] = []
    genotypes: Dict[Tuple[int, str, str], str] = {}

    for line in vcf_content.splitlines():
        if not line:
            continue
        if line.startswith("#"):
            header.append(line)
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4].split(",")[0]
        format_fields = fields[8].split(":")
        sample_fields = fields[9].split(":")
        if "GT" not in format_fields:
            continue
        gt_index = format_fields.index("GT")
        genotypes[(pos, ref, alt)] = sample_fields[gt_index]

    return header, genotypes


def _build_vcf_line(call: RefinedCall) -> str:
    qual = min(999, int(call.confidence * 100))
    info = f"DP={call.depth};AF={call.confidence:.3f}"
    return (
        f"{call.chrom}\t{call.pos}\t.\t{call.ref}\t{call.alt}\t{qual}\tPASS\t"
        f"{info}\tGT\t{call.genotype}"
    )


def merge_refined_calls(
    vcf_content: str,
    refined_calls: List[RefinedCall],
    sample_name: str = SAMPLE_NAME,
) -> str:
    """Replace bcftools records with diffusion-refined genotypes."""
    header, _ = _parse_vcf_records(vcf_content)
    if not header:
        header = [
            "##fileformat=VCFv4.2",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Diffusion confidence">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}",
        ]

    for index, line in enumerate(header):
        if line.startswith("#CHROM"):
            fields = line.split("\t")
            fields[-1] = sample_name
            header[index] = "\t".join(fields)

    body = [_build_vcf_line(call) for call in sorted(refined_calls, key=lambda c: c.pos)]
    output = header + body
    text = "\n".join(output)
    return text + ("\n" if text else "")


def refine_vcf_with_diffusion(
    vcf_content: str,
    bam_path: str,
    ref_fasta: str,
    region: str,
    *,
    steps: Optional[int] = None,
    min_confidence: Optional[float] = None,
) -> str:
    """Run diffusion refinement on an existing bcftools VCF."""
    if os.environ.get("NIOME_DISABLE_DIFFUSION", "").lower() in ("1", "true", "yes"):
        return vcf_content

    _, existing_genotypes = _parse_vcf_records(vcf_content)
    pileups = collect_pileup_sites(bam_path, ref_fasta, region)

    refiner = DiffusionGenotypeRefiner(
        steps=steps or int(os.environ.get("NIOME_DIFFUSION_STEPS", DEFAULT_STEPS)),
        min_confidence=min_confidence
        or float(os.environ.get("NIOME_DIFFUSION_MIN_CONF", DEFAULT_MIN_CONFIDENCE)),
    )
    refined = refiner.refine_calls(pileups, existing_genotypes)

    if not refined:
        return vcf_content

    return merge_refined_calls(vcf_content, refined)
