"""Tests for variant calling helpers."""

from niome_subnet.genomics.variant_caller import (
    VariantCallingError,
    _validate_reference_fasta,
    _vcf_column_names,
    annotate_cftr_variants,
    build_empty_vcf,
    ensure_read_file,
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


def test_ensure_read_file_copies_local_fastq(tmp_path):
    source = tmp_path / "reads_1.fq"
    source.write_text("@read1\nACGT\n", encoding="utf-8")
    dest = tmp_path / "work" / "read_1.fq"

    ensure_read_file(str(source), str(dest))

    assert dest.read_text(encoding="utf-8") == "@read1\nACGT\n"


def test_ensure_read_file_skips_when_dest_exists(tmp_path):
    source = tmp_path / "reads_1.fq"
    source.write_text("@read1\nACGT\n", encoding="utf-8")
    dest = tmp_path / "read_1.fq"
    dest.write_text("cached\n", encoding="utf-8")

    ensure_read_file(str(source), str(dest))

    assert dest.read_text(encoding="utf-8") == "cached\n"


def test_ensure_reference_rejects_vcf(tmp_path):
    vcf_path = tmp_path / "truth.vcf"
    vcf_path.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n", encoding="utf-8")

    try:
        _validate_reference_fasta(str(vcf_path))
        assert False, "expected VariantCallingError"
    except VariantCallingError as exc:
        assert "VCF file" in str(exc)


def test_ensure_reference_replaces_invalid_file(tmp_path, monkeypatch):
    ref_path = tmp_path / "ref.fa"
    ref_path.write_text("not fasta\n", encoding="utf-8")

    def fake_fetch(chrom, start, end, dest):
        with open(dest, "w", encoding="utf-8") as handle:
            handle.write(f">{chrom}\n" + ("A" * (end - start)) + "\n")

    monkeypatch.setattr(
        "niome_subnet.genomics.variant_caller.fetch_reference_region",
        fake_fetch,
    )

    path, offset = ensure_reference(
        "chr7", 117_480_000, 117_670_000, ref_fasta=str(ref_path)
    )
    assert path == str(ref_path)
    assert offset == 117_480_000
    assert ref_path.read_text(encoding="utf-8").startswith(">chr7\n")


def test_ensure_reference_detects_regional_fasta(tmp_path):
    ref_path = tmp_path / "ref.fa"
    ref_path.write_text(">chr7\n" + ("A" * 190_000) + "\n", encoding="utf-8")

    _, offset = ensure_reference(
        "chr7", 117_480_000, 117_670_000, ref_fasta=str(ref_path)
    )
    assert offset == 117_480_000


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
