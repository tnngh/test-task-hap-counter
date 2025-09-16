# The task

## Goal

Produce a program in Python programming language, which given read alignments (in BAM format with haplotype tags) and a set of variants (in phased VCF format) computes support for ALT and REF alleles across reads assigned to individual haplotypes.

**Note that usage of GenAI coding assistants (GitHub Copilot, Claude Code, etc) is both allowed and encouraged**.

## Inputs & assumptions

* **alignments.bam** -- BAM file with read alignments. Coordinate-sorted, indexed (matching .bai index is also provided) alignments BAM, with primary alignments optionally **haplotagged** with`HP`tag (1 or 2).
* **variants.vcf(.gz)** -- variants in phased VCF format (optionally gzipped). To keep things relatively simple, we only ask you to consider bi-allelic mismatch variants (SNVs), ignoring indels, MNPs, *etc*.


## Output
 
Produce a TSV file with rows corresponding to SNVs and the following 6 required columns (feel free to include extra columns if you think they might be useful):

* `chrom` – chromosome name
* `pos`– position on chromosome
* `h<H>_<A>` with `H` in {1,2} and `A` in {'ALT, 'REF'} (4 columns total) — number of primary alignments with HP==`H` supporting allele `A`

# Test data

Test data archive contains read alignments (30x coverage) and variants for chr16:28000000-28500000 region of HG002 genome.

# Development Setup

This project uses Docker for consistent development environment.

## Option 1: VSCode Devcontainer (Recommended)

1. Open in VSCode
2. When prompted, click "Reopen in Container"
3. Wait for container to build and dependencies to install
4. Start developing!

## Option 2: Standalone Docker

```bash
# Build the image
docker build -t haplotype-counter .

# Run interactively
docker run -it --rm -v $(pwd):/workspace haplotype-counter

# Or run with specific command
docker run --rm -v $(pwd):/workspace haplotype-counter haplotype-counter --help
```

## Available Commands

```bash
# Install dependencies (if not already installed)
poetry install

# Run tests
poetry run pytest

# Run tests with coverage
poetry run pytest --cov=haplotype_counter --cov-report=html

# Format code
poetry run black .

# Lint code
poetry run ruff check .

# Type check
poetry run mypy src/

# Run the CLI (once implemented)
poetry run haplotype-counter test_data/*.bam test_data/*.vcf.gz -o results.tsv

# Or activate shell and run directly
poetry shell
haplotype-counter test_data/*.bam test_data/*.vcf.gz -o results.tsv
```

## Examining Test Data

The Docker image includes samtools, bcftools, and other bioinformatics tools:

```bash
# Examine BAM file
samtools view -H test_data/*.bam | head -20
samtools view test_data/*.bam | head -10
samtools view test_data/*.bam | cut -f12- | grep -o 'HP:i:[12]' | sort | uniq -c

# Examine VCF file  
zcat test_data/*.vcf.gz | head -50 | grep "^#"
zcat test_data/*.vcf.gz | grep -v "^#" | head -10

# Check file integrity
samtools quickcheck test_data/*.bam
tabix -l test_data/*.vcf.gz
```

## Project Structure

```
src/haplotype_counter/    # Main package
├── __init__.py
├── cli.py               # Command-line interface
├── core.py             # Core processing logic (to be implemented)
└── utils.py            # Utility functions (to be implemented)

tests/                   # Test suite
├── __init__.py
├── test_cli.py         # CLI tests
└── test_core.py        # Core logic tests (to be implemented)

```