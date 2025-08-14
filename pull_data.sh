# Pick a compact test window
REGION="chr16:28000000-28500000"

VCF_URL="https://42basepairs.com/download/s3/ont-open-data/giab_2023.05/analysis/variant_calling/hg002_sup_60x/hg002.wf_snp.vcf.gz"
BAM_URL="https://42basepairs.com/download/s3/ont-open-data/giab_2023.05/analysis/variant_calling/hg002_sup_60x/hg002.haplotagged.bam"

rm -rf test_data/
mkdir -p test_data/

# Subset the VCF to the region (keep GT/PS; drop heavy INFO/FORMAT fields)
bcftools view -r "$REGION" -Ou "$VCF_URL" \
| bcftools annotate -x INFO,FORMAT/DP,FORMAT/AF,FORMAT/PL,FORMAT/AD \
  -Oz -o test_data/giab_2023.05.hg002.wf_snp.${REGION//[:\-]/_}.vcf.gz

# Index the regional VCF
bcftools index -t -f test_data/giab_2023.05.hg002.wf_snp.${REGION//[:\-]/_}.vcf.gz

# Write BAM + index in one go (samtools '##idx##' syntax)
# Use -s 42.5 to subsample the BAM to 30x coverage
samtools view -@ 16 -b "$BAM_URL" "$REGION" -s 42.5 --write-index \
  -o test_data/giab_2023.05.hg002.haplotagged.${REGION//[:\-]/_}.30x.bam##idx##test_data/giab_2023.05.hg002.haplotagged.${REGION//[:\-]/_}.30x.bam.bai

# Remove indexes from the original files
rm -f hg002.haplotagged.bam.bai hg002.wf_snp.vcf.gz.tbi