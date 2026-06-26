"""Extract per-position pileup features from aligned reads."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

import pysam

_BASES = ("A", "C", "G", "T")


@dataclass
class PileupSite:
    chrom: str
    pos: int
    ref: str
    depth: int
    base_counts: Dict[str, int] = field(default_factory=dict)
    alt_allele: Optional[str] = None
    alt_fraction: float = 0.0

    def alt_support(self) -> int:
        if self.alt_allele is None:
            return 0
        return self.base_counts.get(self.alt_allele.upper(), 0)


def _parse_region(region: str) -> Tuple[str, int, int]:
    chrom, bounds = region.split(":")
    start_s, end_s = bounds.split("-")
    return chrom, int(start_s), int(end_s)


def _load_reference_bases(ref_fasta: str, chrom: str, start: int, end: int) -> Dict[int, str]:
    bases: Dict[int, str] = {}
    with pysam.FastaFile(ref_fasta) as fasta:
        if chrom not in fasta.references:
            raise ValueError(f"Reference missing contig {chrom}")
        for pos in range(start, end + 1):
            bases[pos] = fasta.fetch(chrom, pos - 1, pos).upper()
    return bases


def _dominant_alt(ref_base: str, base_counts: Dict[str, int]) -> Tuple[Optional[str], float]:
    ref = ref_base.upper()
    depth = sum(base_counts.values())
    if depth == 0:
        return None, 0.0

    best_alt = None
    best_count = 0
    for base, count in base_counts.items():
        if base == ref:
            continue
        if count > best_count:
            best_alt = base
            best_count = count

    if best_alt is None:
        return None, 0.0
    return best_alt, best_count / depth


def iter_pileup_sites(
    bam_path: str,
    ref_fasta: str,
    region: str,
    *,
    min_depth: int = 6,
    max_depth: int = 500,
) -> Iterator[PileupSite]:
    """Yield pileup features for positions with sufficient coverage."""
    chrom, start, end = _parse_region(region)
    ref_bases = _load_reference_bases(ref_fasta, chrom, start, end)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for column in bam.pileup(
            chrom,
            start - 1,
            end,
            truncate=True,
            stepper="all",
            min_base_quality=20,
        ):
            pos = column.reference_pos + 1
            ref_base = ref_bases.get(pos, "N")
            if ref_base not in _BASES:
                continue

            counts: Dict[str, int] = {base: 0 for base in _BASES}
            for pileup_read in column.pileups:
                if pileup_read.is_refskip or pileup_read.is_del:
                    continue
                base = pileup_read.query_sequence[pileup_read.query_position].upper()
                if base in counts:
                    counts[base] += 1

            depth = sum(counts.values())
            if depth < min_depth or depth > max_depth:
                continue

            alt, alt_fraction = _dominant_alt(ref_base, counts)
            yield PileupSite(
                chrom=chrom,
                pos=pos,
                ref=ref_base,
                depth=depth,
                base_counts=counts,
                alt_allele=alt,
                alt_fraction=alt_fraction,
            )


def collect_pileup_sites(
    bam_path: str,
    ref_fasta: str,
    region: str,
    **kwargs,
) -> List[PileupSite]:
    return list(iter_pileup_sites(bam_path, ref_fasta, region, **kwargs))
