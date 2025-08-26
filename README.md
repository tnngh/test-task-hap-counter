# The task

## Goal

Given:

* a **haplotagged BAM** (primary alignments carry `HP` tag),
* a **VCF/VCF.gz** of variants,

produce a DataFrame with one row per SNV and four required integer columns:

* `chrom`, `pos`, `h1_ref`, `h1_alt`, `h2_ref`, `h2_alt`

Optionally, you may include extra columns like if convenient and find useful.

## Inputs & assumptions

* BAM is coordinate-sorted, indexed, and **haplotagged** via `HP` tag (1 or 2).
* VCF is phased and may contain mixed variant types; **process only rows where `len(REF)==1` and `len(ALT)==1`** (SNVs). You can ignore indels, MNPs, *etc*.

## Output schema

Each row corresponds to one **input VCF SNP** that passed filters. Required columns:

* `chrom` – chromosome
* `pos` – position on chromosome
* `h1_ref` — number of HP=1 reads supporting REF
* `h1_alt` — number of HP=1 reads supporting ALT
* `h2_ref` — number of HP=2 reads supporting REF
* `h2_alt` — number of HP=2 reads supporting ALT

# Test data

Test data contains a haplotagged BAM file (30x coverage) and corresponding VCF with variants for region chr16:28000000-28500000 of HG002.
You can run `pull_data.sh` to download the data from GIAB HG002 dataset or pull it from this repository.