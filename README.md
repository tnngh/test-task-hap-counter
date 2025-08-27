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
