"""Tests for diffusion genotype refiner."""

import numpy as np

from niome_subnet.genomics.ml.diffusion_refiner import (
    DiffusionGenotypeRefiner,
    RefinedCall,
    merge_refined_calls,
)
from niome_subnet.genomics.ml.pileup import PileupSite


def _site(alt_fraction: float, depth: int = 20, ref: str = "C", alt: str = "A") -> PileupSite:
    return PileupSite(
        chrom="chr7",
        pos=117509084,
        ref=ref,
        depth=depth,
        base_counts={ref: int(depth * (1 - alt_fraction)), alt: int(depth * alt_fraction)},
        alt_allele=alt,
        alt_fraction=alt_fraction,
    )


def test_denoise_prefers_het_at_half_alt_fraction():
    refiner = DiffusionGenotypeRefiner(steps=20)
    posterior = refiner.denoise(np.array([0.34, 0.33, 0.33]), alt_fraction=0.5, depth=20)
    assert posterior[1] == max(posterior)


def test_refine_site_filters_low_confidence():
    refiner = DiffusionGenotypeRefiner(min_confidence=0.9, min_alt_fraction=0.05)
    call = refiner.refine_site(_site(0.08, depth=8))
    assert call is None


def test_refine_site_calls_het_variant():
    refiner = DiffusionGenotypeRefiner(min_confidence=0.5, min_alt_fraction=0.15)
    call = refiner.refine_site(_site(0.45, depth=24), prior_gt="0/1")
    assert call is not None
    assert call.genotype in ("0/1", "1/1")
    assert call.pos == 117509084


def test_merge_refined_calls_preserves_header():
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    merged = merge_refined_calls(
        vcf,
        [
            RefinedCall(
                chrom="chr7",
                pos=117509084,
                ref="C",
                alt="A",
                genotype="0/1",
                confidence=0.82,
                depth=22,
            )
        ],
    )
    assert "117509084" in merged
    assert "0/1" in merged
    assert merged.startswith("##fileformat=VCFv4.2")
