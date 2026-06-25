"""Tests for local ground-truth comparison helpers."""

from niome_subnet.genomics.compare import (
    _hgvs_matches_variant,
    map_vcf_to_cftr2_annotations,
)


TRUTH_ANNOTATIONS = {
    "4064589": {
        "hgvs": "NC_000007.14:g.117509084C>A",
        "clinical_significance": "Likely pathogenic",
        "drug_response": {
            "ivacaftor": "non_responsive",
            "tezacaftor_ivacaftor": "non_responsive",
            "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
            "lumacaftor_ivacaftor": "non_responsive",
        },
    },
    "1772239": {
        "hgvs": "NC_000007.14:g.117504341_117504344del",
        "clinical_significance": "Pathogenic",
        "drug_response": {
            "ivacaftor": "non_responsive",
            "tezacaftor_ivacaftor": "non_responsive",
            "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
            "lumacaftor_ivacaftor": "non_responsive",
        },
    },
}


TRUTH_VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr7\t117509084\t.\tC\tA\t.\tPASS\t.\tGT\t0/1
chr7\t117504339\t.\tACAAT\tA\t.\tPASS\t.\tGT\t1/1
"""


def test_hgvs_matches_snv():
    assert _hgvs_matches_variant("NC_000007.14:g.117509084C>A", 117509084, "C", "A")


def test_hgvs_matches_deletion_range():
    assert _hgvs_matches_variant(
        "NC_000007.14:g.117504341_117504344del", 117504339, "ACAAT", "A"
    )


def test_map_vcf_to_cftr2_annotations(tmp_path):
    truth_vcf = tmp_path / "truth.vcf"
    truth_vcf.write_text(TRUTH_VCF, encoding="utf-8")
    annotations_path = tmp_path / "annotations.json"
    annotations_path.write_text(
        __import__("json").dumps(TRUTH_ANNOTATIONS),
        encoding="utf-8",
    )

    miner_vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr7\t117509084\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\n"
    )
    mapped = map_vcf_to_cftr2_annotations(
        miner_vcf,
        truth_vcf_path=str(truth_vcf),
        truth_annotations_path=str(annotations_path),
    )

    assert set(mapped) == {"4064589"}
    assert mapped["4064589"]["hgvs"] == "NC_000007.14:g.117509084C>A"
    assert mapped["4064589"]["clinical_significance"] == "unknown"
