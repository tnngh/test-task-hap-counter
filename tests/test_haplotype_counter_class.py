"""Tests for the HaplotypeCounter class."""
from pathlib import Path
from re import M
from unittest.mock import Mock, patch

import pytest

from haplotype_counter.haplotype_counter import SNV, HaplotypeCount, HaplotypeCounter


def test_haplotype_counter_initialization() -> None:
    """Test HaplotypeCounter initialization with valid parameters."""
    bam_file = Path("test.bam")

    with patch.object(Path, 'exists', return_value=True):
        counter = HaplotypeCounter(
            bam_file=bam_file,
            min_base_quality=20,
            min_mapping_quality=30,
            max_depth=500,
            stepper='nofilter'
        )

        assert counter._bam_file == bam_file
        assert counter._min_base_quality == 20
        assert counter._min_mapping_quality == 30
        assert counter._max_depth == 500
        assert counter._stepper == 'nofilter'
        assert counter._bam_handle is None


def test_haplotype_counter_bam_file_not_found() -> None:
    """Test HaplotypeCounter raises FileNotFoundError when BAM file doesn't exist."""
    bam_file = Path("nonexistent.bam")

    with patch.object(Path, 'exists', return_value=False):
        with pytest.raises(FileNotFoundError, match="BAM file not found"):
            HaplotypeCounter(bam_file=bam_file)


def test_haplotype_counter_bai_file_not_found() -> None:
    """Test HaplotypeCounter raises FileNotFoundError when BAI file doesn't exist."""
    # Create a temporary BAM file that will be deleted after the test
    import tempfile

    with tempfile.NamedTemporaryFile(suffix=".bam") as temp_bam:
        bam_file = Path(temp_bam.name)

        # The BAM file exists (created by tempfile) but the BAI file doesn't
        with pytest.raises(FileNotFoundError, match="BAM index file not found"):
            with HaplotypeCounter(bam_file=bam_file):
                pass


def test_haplotype_counter_context_manager() -> None:
    """Test HaplotypeCounter as context manager."""
    bam_file = Path("test.bam")

    with patch.object(Path, 'exists', return_value=True):
        with patch('pysam.AlignmentFile') as mock_alignment_file:
            mock_bam = Mock()
            mock_alignment_file.return_value = mock_bam

            with HaplotypeCounter(bam_file=bam_file) as counter:
                assert counter._bam_handle == mock_bam
                mock_alignment_file.assert_called_once_with(str(bam_file), "rb")

            # Verify cleanup
            mock_bam.close.assert_called_once()
            assert counter._bam_handle is None


def test_haplotype_counter_bam_property() -> None:
    """Test HaplotypeCounter bam property opens file when needed."""
    bam_file = Path("test.bam")

    with patch.object(Path, 'exists', return_value=True):
        with patch('pysam.AlignmentFile') as mock_alignment_file:
            mock_bam = Mock()
            mock_alignment_file.return_value = mock_bam

            counter: HaplotypeCounter = HaplotypeCounter(bam_file=bam_file)
            bam_handle = counter.bam

            assert bam_handle == mock_bam
            assert counter._bam_handle == mock_bam
            mock_alignment_file.assert_called_once_with(str(bam_file), "rb")


def _create_mock_pileup_reads(
    query_name: str, query_sequence: str, hp_tag: int, query_position: int
) -> tuple[Mock, Mock]:
    """Create mock PileupRead, AlignedSegment tuple that won't be skipped by the HaplotypeCounter."""
    mock_pileupread = Mock()
    mock_pileupread.is_del = False
    mock_pileupread.is_refskip = False
    mock_pileupread.indel = 0

    mock_alignedsegment = Mock()
    mock_alignedsegment.is_secondary = False
    mock_alignedsegment.is_supplementary = False
    mock_alignedsegment.is_duplicate = False
    mock_alignedsegment.is_qcfail = False
    mock_alignedsegment.query_name = query_name
    mock_alignedsegment.query_sequence = query_sequence
    mock_alignedsegment.get_tag.return_value = hp_tag  # HP tag for haplotype 1

    mock_pileupread.alignment = mock_alignedsegment
    mock_pileupread.query_position = query_position  # Position in read matching reference base

    return mock_pileupread, mock_alignedsegment


def test_haplotype_counter_count_single_snv() -> None:
    """Test HaplotypeCounter count_single_snv method."""
    bam_file = Path("test.bam")
    snv = SNV(chrom="chr1", pos_1based=12345, ref="A", alt="T")

    with patch.object(Path, 'exists', return_value=True):
        with patch('pysam.AlignmentFile') as mock_alignment_file:
            # Mock BAM file and pileup
            mock_bam = Mock()
            mock_alignment_file.return_value = mock_bam

            # Mock pileup column
            mock_pileupcolumn = Mock()
            mock_pileupcolumn.reference_pos = 12344  # 0-based position

            # Create pileup reads for different haplotypes and alleles
            # H1 ref read
            h1_ref_pileupread, h1_ref_alignment = _create_mock_pileup_reads(
                query_name="read1", query_sequence="ATCG", hp_tag=1, query_position=0
            )

            # H1 alt read
            h1_alt_pileupread, h1_alt_alignment = _create_mock_pileup_reads(
                query_name="read2", query_sequence="TTCG", hp_tag=1, query_position=0
            )

            # H2 ref read
            h2_ref_pileupread, h2_ref_alignment = _create_mock_pileup_reads(
                query_name="read3", query_sequence="ATCG", hp_tag=2, query_position=0
            )

            h2_alt_pileupread, h2_alt_alignment = _create_mock_pileup_reads(
                query_name="read4", query_sequence="TTCG", hp_tag=2, query_position=0
            )

            # H2 alt read
            h2_alt_pileupread, h2_alt_alignment = _create_mock_pileup_reads(
                query_name="read4", query_sequence="TTCG", hp_tag=2, query_position=0
            )

            # Add all reads to the pileup
            mock_pileupcolumn.pileups = [h1_ref_pileupread, h1_alt_pileupread, h2_ref_pileupread, h2_alt_pileupread]
            mock_bam.pileup.return_value = [mock_pileupcolumn]

            counter: HaplotypeCounter = HaplotypeCounter(bam_file=bam_file)
            result: HaplotypeCount = counter.count_single_snv(snv)

            assert isinstance(result, HaplotypeCount)
            assert result.chrom == "chr1"
            assert result.pos_1based == 12345
            # Check that all counts are non-zero
            assert result.h1_ref == 1
            assert result.h1_alt == 1
            assert result.h2_ref == 1
            assert result.h2_alt == 1


def test_haplotype_counter_process_snvs() -> None:
    """Test HaplotypeCounter process_snvs method returns a list of HaplotypeCount objects.

    There should be a HaplotypeCount object for each SNV.
    """
    bam_file: Path = Path("test.bam")
    snvs: list[SNV] = [
        SNV(chrom="chr1", pos_1based=12345, ref="A", alt="T"),
        SNV(chrom="chr1", pos_1based=12346, ref="G", alt="C")
    ]

    with patch.object(Path, 'exists', return_value=True):
        with patch('pysam.AlignmentFile') as mock_alignment_file:
            # Mock BAM file and pileup
            mock_bam = Mock()
            mock_alignment_file.return_value = mock_bam

            # Mock pileup columns for each SNV
            mock_pileupcolumn1 = Mock()
            mock_pileupcolumn1.reference_pos = 12344  # 0-based position for first SNV

            mock_pileupcolumn2 = Mock()
            mock_pileupcolumn2.reference_pos = 12345  # 0-based position for second SNV

            # Mock pileup reads for first SNV
            mock_pileupread1_h1, mock_alignment1_h1 = _create_mock_pileup_reads(
                query_name="read1", query_sequence="ATCG", hp_tag=1, query_position=0
            )
            mock_pileupread1_h2, mock_alignment1_h2 = _create_mock_pileup_reads(
                query_name="read2", query_sequence="TTCG", hp_tag=2, query_position=0
            )

            mock_pileupcolumn1.pileups = [mock_pileupread1_h1, mock_pileupread1_h2]

            # Mock pileup reads for second SNV
            mock_pileupread2_h1, mock_alignment2_h1 = _create_mock_pileup_reads(
                query_name="read3", query_sequence="CCTA", hp_tag=1, query_position=0
            )
            mock_pileupread2_h2, mock_alignment2_h2 = _create_mock_pileup_reads(
                query_name="read4", query_sequence="GCTA", hp_tag=2, query_position=0
            )

            mock_pileupcolumn2.pileups = [mock_pileupread2_h1, mock_pileupread2_h2]

            # Set up mock pileup to return our columns
            mock_bam.pileup.side_effect = lambda contig, start, **kwargs: \
                [mock_pileupcolumn1] if contig == "chr1" and start == 12344 else [mock_pileupcolumn2]

            counter: HaplotypeCounter = HaplotypeCounter(bam_file=bam_file)
            results: list[HaplotypeCount] = counter.process_snvs(snvs)

            assert len(results) == 2
            assert all(isinstance(result, HaplotypeCount) for result in results)
            assert results[0].chrom == "chr1"
            assert results[0].pos_1based == 12345
            assert results[0].h1_ref == 1  # One haplotype 1 read with ref allele
            assert results[0].h1_alt == 0  # No haplotype 1 reads with alt allele
            assert results[0].h2_ref == 0  # No haplotype 2 reads with ref allele
            assert results[0].h2_alt == 1  # One haplotype 2 read with alt allele

            assert results[1].chrom == "chr1"
            assert results[1].pos_1based == 12346
            assert results[1].h1_ref == 0  # One haplotype 1 read with ref allele
            assert results[1].h1_alt == 1  # No haplotype 1 reads with alt allele
            assert results[1].h2_ref == 1  # No haplotype 2 reads with ref allele
            assert results[1].h2_alt == 0  # One haplotype 2 read with alt allele


def test_haplotype_counter_process_vcf() -> None:
    """Test HaplotypeCounter process_vcf method."""
    bam_file = Path("test.bam")
    vcf_file = Path("test.vcf")

    with patch.object(Path, 'exists', return_value=True):
        with patch('pysam.AlignmentFile') as mock_alignment_file:
            with patch('haplotype_counter.haplotype_counter.parse_vcf_snvs') as mock_parse:
                # Mock VCF parsing
                mock_snvs = [SNV(chrom="chr1", pos_1based=12345, ref="A", alt="T")]
                mock_parse.return_value = mock_snvs

                # Mock BAM file and pileup
                mock_bam = Mock()
                mock_alignment_file.return_value = mock_bam
                mock_bam.pileup.return_value = []

                counter: HaplotypeCounter = HaplotypeCounter(bam_file=bam_file)
                results: list[HaplotypeCount] = counter.process_vcf(vcf_file)

                mock_parse.assert_called_once_with(vcf_file)
                assert len(results) == 1
                assert isinstance(results[0], HaplotypeCount)
                assert results[0].chrom == "chr1"
                assert results[0].pos_1based == 12345
