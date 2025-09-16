"""Core functionality for counting haplotype allele support."""

import logging
from array import array
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pysam

logger: logging.Logger = logging.getLogger(__name__)


@dataclass
class SNV:
    """Represents a single nucleotide variant.

    Args:
        chrom: Chromosome name
        pos: Position on chromosome (1-based)
        ref: Reference allele
        alt: Alternate allele
    """

    chrom: str
    pos_1based: int
    ref: str
    alt: str


@dataclass
class HaplotypeCount:
    """Represents allele counts for a specific SNV across haplotypes.

    Assumes:
     - organism has max 2 haplotypes

    Args:
        chrom: Chromosome name
        pos: Position on chromosome (1-based)
        h1_ref: Number of reads supporting the REF allele for haplotype 1
        h1_alt: Number of reads supporting the ALT allele for haplotype 1
        h2_ref: Number of reads supporting the REF allele for haplotype 2
        h2_alt: Number of reads supporting the ALT allele for haplotype 2
    """

    chrom: str
    pos_1based: int
    h1_ref: int = 0
    h1_alt: int = 0
    h2_ref: int = 0
    h2_alt: int = 0


class HaplotypeCounter:
    """Haplotype allele counting.

    This class provides efficient processing of haplotype allele counts by
    maintaining open file handles and configuration state across multiple
    SNV processing calls.

    Args:
        bam_file: Path to indexed BAM file with haplotagged reads
        min_base_quality: Minimum base quality threshold (default: 0)
        min_mapping_quality: Minimum mapping quality threshold (default: 0)
        max_depth: Maximum depth limit for pileup (default: 1000)
        stepper: Pileup stepper mode (default: 'nofilter')
    """

    def __init__(
        self,
        bam_file: Path,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
        max_depth: int = 1000,
        stepper: str = 'nofilter'
    ) -> None:
        """Initialize the haplotype counter with configuration."""
        self._bam_file: Path = bam_file
        self._min_base_quality: int = min_base_quality
        self._min_mapping_quality = min_mapping_quality
        self._max_depth: int = max_depth
        self._stepper: str = stepper
        self._bam_handle: pysam.AlignmentFile | None = None

        # Validate BAM file exists
        if not self._bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {self._bam_file}")

        # Validate BAM index exists
        bai_file = Path(str(self._bam_file) + ".bai")
        if not bai_file.exists():
            raise FileNotFoundError(f"BAM index file not found: {bai_file}")

    def __enter__(self) -> "HaplotypeCounter":
        """Context manager entry - open BAM file."""
        self._bam_handle = pysam.AlignmentFile(str(self._bam_file), "rb")
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Context manager exit - close BAM file."""
        if self._bam_handle is not None:
            self._bam_handle.close()
            self._bam_handle = None

    @property
    def bam(self) -> pysam.AlignmentFile:
        """Get the BAM file handle, opening it if necessary."""
        if self._bam_handle is None:
            self._bam_handle = pysam.AlignmentFile(str(self._bam_file), "rb")
        return self._bam_handle

    def count_single_snv(self, snv: SNV) -> HaplotypeCount:
        """Count allele support for a single SNV.

        Args:
            snv: SNV to analyze

        Returns:
            HaplotypeCount object with counts for each haplotype and allele

        Raises:
            ValueError: If haplotype tag is invalid
        """
        counts = HaplotypeCount(chrom=snv.chrom, pos_1based=snv.pos_1based)

        for pileupcolumn in self.bam.pileup(
            contig=snv.chrom,
            start=snv.pos_1based - 1,  # pysam uses 0-based position
            end=snv.pos_1based,
            stepper=self._stepper,
            min_base_quality=self._min_base_quality,
            min_mapping_quality=self._min_mapping_quality,
            max_depth=self._max_depth
        ):
            if pileupcolumn.reference_pos == snv.pos_1based - 1:  # pysam uses 0-based position
                for pileupread in pileupcolumn.pileups:
                    if not is_primary_align_qc_pass_mismatch(pileupread):
                        continue

                    read: pysam.AlignedSegment = pileupread.alignment

                    query_position: int | None = pileupread.query_position
                    if query_position is None or read.query_sequence is None:
                        logger.warning(
                            f"Skipping read {read.query_name} because it has no query position or query sequence")
                        continue

                    base: str = read.query_sequence[query_position]

                    # Get haplotype tag
                    try:
                        haplotype: str | int | float | array = read.get_tag(tag="HP")
                    except KeyError:
                        continue

                    if base is None or haplotype is None or not isinstance(haplotype, int):
                        logger.warning(
                            f"Skipping read {read.query_name} because it has "
                            f"invalid base or haplotype tag at ref position {snv.pos_1based}")
                        continue

                    # Count allele support
                    if haplotype not in [1, 2]:
                        raise ValueError(
                            f"Invalid haplotype tag: {haplotype} for read {read.query_name}.  "
                            "Expected 1 or 2.")

                    if base.upper() == snv.ref.upper():
                        if haplotype == 1:
                            counts.h1_ref += 1
                        elif haplotype == 2:
                            counts.h2_ref += 1
                    elif base.upper() == snv.alt.upper():
                        if haplotype == 1:
                            counts.h1_alt += 1
                        elif haplotype == 2:
                            counts.h2_alt += 1

        return counts

    def process_snvs(self, snvs: list[SNV]) -> list[HaplotypeCount]:
        """Process multiple SNVs.

        Args:
            snvs: List of SNVs to process

        Returns:
            List of HaplotypeCount objects, one for each SNV
        """
        results: list[HaplotypeCount] = []
        for snv in snvs:
            count = self.count_single_snv(snv)
            results.append(count)
        return results

    def process_vcf(self, vcf_file: Path) -> list[HaplotypeCount]:
        """Parse VCF file and process all SNVs.

        Args:
            vcf_file: Path to VCF file with variants

        Returns:
            List of HaplotypeCount objects, one for each SNV

        Raises:
            FileNotFoundError: If VCF file doesn't exist
            ValueError: If VCF file is malformed
        """
        snvs = parse_vcf_snvs(vcf_file)
        return self.process_snvs(snvs)


def parse_vcf_snvs(vcf_file: Path) -> list[SNV]:
    """Parse VCF file and extract bi-allelic SNVs.

    Ideally, we would make this an iterator since VCFs can be large,
    but let's keep it simple for now.

    Args:
        vcf_file: Path to VCF file (can be gzipped)

    Returns:
        List of SNV objects representing bi-allelic single nucleotide variants

    Raises:
        FileNotFoundError: If VCF file doesn't exist
        ValueError: If VCF file is malformed
    """
    if not vcf_file.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_file}")

    snvs: list[SNV] = []

    try:
        with pysam.VariantFile(str(vcf_file)) as vcf:
            for record in vcf:
                # Skip non-bi-allelic variants
                if record.alts is None or len(record.alts) != 1:
                    continue

                # Skip indels - only keep SNVs (single base changes)
                ref: str | None = record.ref
                alt: str | None = record.alts[0]

                if ref is not None and alt is not None and len(ref) == 1 and len(alt) == 1:
                    snvs.append(
                        SNV(
                            chrom=record.chrom,
                            pos_1based=record.pos,  # pysam uses 0-based position
                            ref=ref,
                            alt=alt,
                        )
                    )
    except Exception as e:
        raise ValueError(
            f"Error parsing VCF file {vcf_file}: {str(e)}"
        ) from e

    return snvs


def is_primary_align_qc_pass_mismatch(pileupread: pysam.PileupRead) -> bool:
    """Check if an alignment is a primary alignment mismatch that passes QC.

    Args:
        read: pysam PileupRead object

    Returns:
        True if:
        - the alignment is a primary alignment AND
        - the alignment is non-chimeric AND
        - the alignment is not an optical duplicate AND
        - the alignment is a mismatch (SNV) AND
        - the alignment has a QC pass
        False otherwise
    """
    # Skip deletions, reference skips at this position, and insertions
    if pileupread.is_del or pileupread.is_refskip or pileupread.indel > 0:
        return False

    read: pysam.AlignedSegment = pileupread.alignment

    # Ignore any non-primary alignments, supplementary alignments (eg chimeric alignments),
    # duplicate alignments, and alignments that failed QC.
    if read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
        return False

    return True


def export_to_tsv(counts: list[HaplotypeCount], output_file: Path) -> None:
    """Export haplotype counts to TSV file.

    Args:
        counts: List of HaplotypeCount objects
        output_file: Path to output TSV file

    Raises:
        ValueError: If no counts provided
    """
    if not counts:
        raise ValueError("No counts provided for export")

    # Convert to DataFrame
    data: list[dict[str, str | int]] = []
    for count in counts:
        data.append(
            {
                "chrom": count.chrom,
                "pos": count.pos_1based,
                "h1_REF": count.h1_ref,
                "h1_ALT": count.h1_alt,
                "h2_REF": count.h2_ref,
                "h2_ALT": count.h2_alt,
            }
        )

    df: pd.DataFrame = pd.DataFrame(data)

    # Export to TSV
    df.to_csv(output_file, sep="\t", index=False)
