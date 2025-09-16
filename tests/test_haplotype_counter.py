"""Tests for core haplotype counting functionality."""
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from haplotype_counter.haplotype_counter import (
    SNV,
    HaplotypeCount,
    HaplotypeStats,
    count_allele_support_pileup,
    export_to_tsv,
    get_summary_stats,
    is_primary_align_qc_pass_mismatch,
    parse_vcf_snvs,
    process_all_snvs,
)


def test_snv_creation() -> None:
    """Test SNV dataclass creation."""
    snv: SNV = SNV(chrom="chr1", pos_1based=12345, ref="A", alt="T")
    assert snv.chrom == "chr1"
    assert snv.pos_1based == 12345
    assert snv.ref == "A"
    assert snv.alt == "T"


def test_haplotype_count_creation() -> None:
    """Test HaplotypeCount dataclass creation."""
    count: HaplotypeCount = HaplotypeCount(
        chrom="chr1",
        pos_1based=12345,
        h1_ref=5,
        h1_alt=3,
        h2_ref=7,
        h2_alt=2
    )
    assert count.chrom == "chr1"
    assert count.pos_1based == 12345
    assert count.h1_ref == 5
    assert count.h1_alt == 3
    assert count.h2_ref == 7
    assert count.h2_alt == 2


def test_haplotype_count_defaults() -> None:
    """Test HaplotypeCount with default values."""
    count: HaplotypeCount = HaplotypeCount(chrom="chr1", pos_1based=12345)
    assert count.h1_ref == 0
    assert count.h1_alt == 0
    assert count.h2_ref == 0
    assert count.h2_alt == 0


def test_parse_vcf_snvs_file_not_found() -> None:
    """Test parse_vcf_snvs raises FileNotFoundError for missing file."""
    with pytest.raises(FileNotFoundError, match="VCF file not found"):
        parse_vcf_snvs(Path("nonexistent.vcf"))


@patch("pathlib.Path.exists", return_value=True)
@patch("pysam.VariantFile")
def test_parse_vcf_snvs_filters_correctly(mock_variant_file: Mock, mock_exists: Mock) -> None:
    """Test that parse_vcf_snvs filters for bi-allelic SNVs only."""
    # Mock VCF records
    mock_record1 = Mock()
    mock_record1.chrom = "chr1"
    mock_record1.pos = 100
    mock_record1.ref = "A"
    mock_record1.alts = ("T",)  # Bi-allelic SNV

    mock_record2 = Mock()
    mock_record2.alts = ("T", "G")  # Multi-allelic - should be skipped

    mock_record3 = Mock()
    mock_record3.chrom = "chr1"
    mock_record3.pos = 200
    mock_record3.ref = "AT"  # Indel - should be skipped
    mock_record3.alts = ("A",)

    mock_record4 = Mock()
    mock_record4.chrom = "chr1"
    mock_record4.pos = 300
    mock_record4.ref = "G"
    mock_record4.alts = ("C",)  # Bi-allelic SNV

    mock_vcf = Mock()
    mock_vcf.__iter__ = Mock(
        return_value=iter([
            mock_record1, mock_record2, mock_record3, mock_record4
        ])
    )
    mock_variant_file.return_value.__enter__ = Mock(return_value=mock_vcf)
    mock_variant_file.return_value.__exit__ = Mock(return_value=None)

    snvs: list[SNV] = parse_vcf_snvs(Path("test.vcf"))

    # Should only return the two bi-allelic SNVs
    assert len(snvs) == 2
    assert snvs[0].chrom == "chr1"
    assert snvs[0].pos_1based == 100
    assert snvs[0].ref == "A"
    assert snvs[0].alt == "T"
    assert snvs[1].pos_1based == 300
    assert snvs[1].ref == "G"
    assert snvs[1].alt == "C"


def test_get_summary_stats_empty_list() -> None:
    """Test get_summary_stats with empty list."""
    stats: HaplotypeStats = get_summary_stats([])
    assert stats.total_snvs == 0
    assert stats.total_h1_reads == 0
    assert stats.total_h2_reads == 0
    assert stats.snvs_with_h1_support == 0
    assert stats.snvs_with_h2_support == 0
    assert stats.snvs_with_both_haplotypes == 0
    assert stats.avg_h1_reads_per_snv == 0
    assert stats.avg_h2_reads_per_snv == 0


def test_get_summary_stats_with_data() -> None:
    """Test get_summary_stats with sample data."""
    counts: list[HaplotypeCount] = [
        HaplotypeCount(
            chrom="chr1", pos_1based=100, h1_ref=5, h1_alt=3, h2_ref=2, h2_alt=4
        ),
        HaplotypeCount(
            chrom="chr1", pos_1based=200, h1_ref=0, h1_alt=0, h2_ref=8, h2_alt=1
        ),
        HaplotypeCount(
            chrom="chr1", pos_1based=300, h1_ref=6, h1_alt=2, h2_ref=0, h2_alt=0
        ),
    ]

    stats: HaplotypeStats = get_summary_stats(counts)

    assert stats.total_snvs == 3
    assert stats.total_h1_reads == 16  # 5+3+0+0+6+2
    assert stats.total_h2_reads == 15  # 2+4+8+1+0+0
    assert stats.snvs_with_h1_support == 2  # positions 100 and 300
    assert stats.snvs_with_h2_support == 2  # positions 100 and 200
    assert stats.snvs_with_both_haplotypes == 1  # only position 100
    assert stats.avg_h1_reads_per_snv == 5  # 16/3 rounded to int
    assert stats.avg_h2_reads_per_snv == 5  # 15/3 rounded to int


def test_export_to_tsv_empty_list() -> None:
    """Test export_to_tsv raises ValueError with empty list."""
    with pytest.raises(ValueError, match="No counts provided for export"):
        export_to_tsv([], Path("output.tsv"))


def test_export_to_tsv_valid_output() -> None:
    """Test export_to_tsv creates valid output file"""
    counts: list[HaplotypeCount] = [
        HaplotypeCount(
            chrom="chr1", pos_1based=100, h1_ref=5, h1_alt=3, h2_ref=2, h2_alt=4
        ),
        HaplotypeCount(
            chrom="chr1", pos_1based=200, h1_ref=0, h1_alt=6, h2_ref=1, h2_alt=7
        ),
    ]

    # Use a temporary file that will be automatically deleted
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=True) as temp_file:
        temp_path = Path(temp_file.name)
        export_to_tsv(counts, temp_path)

        # Verify file was created
        assert temp_path.exists(), "Output TSV file should be created"

        # Verify file contents
        with open(temp_path, "r") as f:
            lines = f.readlines()
            assert lines[0] == "chrom\tpos\th1_REF\th1_ALT\th2_REF\th2_ALT\n"
            assert lines[1] == "chr1\t100\t5\t3\t2\t4\n"
            assert lines[2] == "chr1\t200\t0\t6\t1\t7\n"


def test_parse_vcf_snvs_real_file() -> None:
    """Test parse_vcf_snvs with a real VCF file."""
    vcf_file = Path("tests/data/test_variants.vcf")
    snvs = parse_vcf_snvs(vcf_file)

    # Should only return bi-allelic SNVs (positions 28000000, 28000001, 28000002)
    # Skip indel (28000003) and multi-allelic (28000004)
    assert len(snvs) == 3

    # Check first SNV
    assert snvs[0].chrom == "16"
    assert snvs[0].pos_1based == 28000000
    assert snvs[0].ref == "A"
    assert snvs[0].alt == "T"

    # Check second SNV
    assert snvs[1].chrom == "16"
    assert snvs[1].pos_1based == 28000001
    assert snvs[1].ref == "G"
    assert snvs[1].alt == "C"

    # Check third SNV
    assert snvs[2].chrom == "16"
    assert snvs[2].pos_1based == 28000002
    assert snvs[2].ref == "T"
    assert snvs[2].alt == "A"


def test_parse_vcf_snvs_empty_file() -> None:
    """Test parse_vcf_snvs with an empty VCF file."""
    # Create an empty VCF file
    empty_vcf = Path("tests/data/empty.vcf")
    with open(empty_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    snvs = parse_vcf_snvs(empty_vcf)
    assert len(snvs) == 0


def test_parse_vcf_snvs_malformed_file() -> None:
    """Test parse_vcf_snvs with a malformed VCF file."""
    # Create a malformed VCF file
    malformed_vcf = Path("tests/data/malformed.vcf")
    with open(malformed_vcf, "w") as f:
        f.write("This is not a valid VCF file\n")
        f.write("Invalid content\n")

    with pytest.raises(ValueError, match="Error parsing VCF file"):
        parse_vcf_snvs(malformed_vcf)


# Tests for is_primary_align_qc_pass_mismatch function
def test_is_primary_align_qc_pass_mismatch_deletion() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for deletions."""
    # Mock pileupread with deletion
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = True
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_refskip() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for reference skips."""
    # Mock pileupread with reference skip
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = True
    mock_pileupread.indel = 0

    # Mock alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_insertion() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for insertions."""
    # Mock pileupread with insertion
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 2  # Positive indel indicates insertion

    # Mock alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_secondary_alignment() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for secondary alignments."""
    # Mock pileupread with valid mismatch
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock secondary alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = True
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_supplementary_alignment() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for supplementary alignments."""
    # Mock pileupread with valid mismatch
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock supplementary alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = True
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_duplicate_alignment() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for duplicate alignments."""
    # Mock pileupread with valid mismatch
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock duplicate alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = True
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_qcfail_alignment() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False for QC failed alignments."""
    # Mock pileupread with valid mismatch
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock QC failed alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = True
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_valid_mismatch() -> None:
    """Test is_primary_align_qc_pass_mismatch returns True for valid mismatch."""
    # Mock pileupread with valid mismatch
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock valid primary alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is True


def test_is_primary_align_qc_pass_mismatch_multiple_failures() -> None:
    """Test is_primary_align_qc_pass_mismatch returns False when multiple conditions fail."""
    # Mock pileupread with deletion AND secondary alignment
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = True
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    # Mock secondary alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = True
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is False


def test_is_primary_align_qc_pass_mismatch_negative_indel() -> None:
    """Test is_primary_align_qc_pass_mismatch returns True for negative indel (deletion in read)."""
    # Mock pileupread with negative indel (deletion in read, not at position)
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = -1  # Negative indel (deletion in read)

    # Mock valid primary alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is True


def test_is_primary_align_qc_pass_mismatch_zero_indel() -> None:
    """Test is_primary_align_qc_pass_mismatch returns True for zero indel (SNV)."""
    # Mock pileupread with zero indel (SNV)
    mock_pileupread: Mock = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0  # Zero indel (SNV)

    # Mock valid primary alignment
    mock_alignment: Mock = Mock()
    mock_alignment.is_secondary = False
    mock_alignment.is_supplementary = False
    mock_alignment.is_duplicate = False
    mock_alignment.is_qcfail = False
    mock_pileupread.alignment = mock_alignment

    result: bool = is_primary_align_qc_pass_mismatch(mock_pileupread)
    assert result is True
