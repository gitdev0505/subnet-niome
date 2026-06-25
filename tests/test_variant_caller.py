"""Tests for variant calling helpers."""

from niome_subnet.genomics.variant_caller import (
    _vcf_column_names,
    annotate_cftr_variants,
    build_empty_vcf,
    ensure_reference,
    normalize_vcf_sample_column,
    parse_region,
)


def test_parse_region():
    chrom, start, end = parse_region("chr7:117480000-117670000")
    assert chrom == "chr7"
    assert start == 117480000
    assert end == 117670000


def test_build_empty_vcf_has_required_columns():
    vcf = build_empty_vcf("chr7:117480000-117670000")
    header = [line for line in vcf.splitlines() if line.startswith("#CHROM")][0]
    columns = header.lstrip("#").split("\t")
    assert columns == [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE",
    ]
    assert not any(line and not line.startswith("#") for line in vcf.splitlines())


def test_annotate_cftr_variants_empty_for_header_only_vcf():
    vcf = build_empty_vcf("chr7:117480000-117670000")
    assert annotate_cftr_variants(vcf) == {}


def test_annotate_cftr_variants_from_vcf_record():
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr7\t117509084\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\n"
    )
    annotations = annotate_cftr_variants(vcf)
    assert len(annotations) == 1
    entry = next(iter(annotations.values()))
    assert entry["hgvs"] == "NC_000007.14:g.117509084C>A"
    assert entry["clinical_significance"] == "unknown"


def test_vcf_column_names_strips_hash_from_chrom():
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
    assert _vcf_column_names(header) == {
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        "SAMPLE",
    }


def test_ensure_reference_detects_regional_fasta(tmp_path):
    ref_path = tmp_path / "ref.fa"
    ref_path.write_text(">chr7\n" + ("A" * 190_000) + "\n", encoding="utf-8")

    _, offset = ensure_reference(
        "chr7", 117_480_000, 117_670_000, ref_fasta=str(ref_path)
    )
    assert offset == 117_479_999


def test_ensure_reference_full_genome_has_no_offset(tmp_path):
    ref_path = tmp_path / "ref.fa"
    ref_path.write_text(">chr7\n" + ("A" * 1_000_000) + "\n", encoding="utf-8")

    _, offset = ensure_reference(
        "chr7", 117_480_000, 117_670_000, ref_fasta=str(ref_path)
    )
    assert offset == 0


def test_normalize_vcf_sample_column_renames_bcftools_default():
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\taligned.bam\n"
    )
    normalized = normalize_vcf_sample_column(vcf)
    header = [line for line in normalized.splitlines() if line.startswith("#CHROM")][0]
    assert header.endswith("\tSAMPLE")
