---
title: 'How-to: call variants from DNA-seq data '
---

https://software.broadinstitute.org/gatk/documentation/quickstart
https://software.broadinstitute.org/gatk/download/
https://software.broadinstitute.org/gatk/best-practices/
# Install GATK4 
Download the precompiled jar binaries. Set paths. 

# Running GATK
```
gatk HaplotypeCaller \
   -R genome.fa \
   -I sample.sort.dedup.realign.bam \
   -O sample.raw_variants.vcf 
```
Or
```
gatk HaplotypeCaller \
   -R genome.fa \
   -I sample.sort.dedup.realign.bam \
   -O sample.raw_variants_gvcf.vcf 
   -ERC GVCF
```


# Split variants into SNPs and indels 
```
gatk SelectVariants \
-R genome.fa  \
-V sample.raw_variants.vcf \
-O raw_indels.vcf \
--select-type-to-include INDEL
 
gatk SelectVariants \
-R genome.fa  \
-V sample.raw_variants.vcf \
-O raw_indels.vcf \
--select-type-to-include SNP
```
 
# Filter variants 
```
gatk VariantFiltration \
-R genome.fa  \
-V raw_snps.vcf \
-O filtered_snps.vcf \
--filter-name "basic_snp_filter1" \
--filter-expression "QD < 2.0" \
--filter-name "basic_snp_filter2" \
--filter-expression "FS > 60.0" \
--filter-name "basic_snp_filter3" \
--filter-expression "MQ < 40.0" \
--filter-name "basic_snp_filter4" \
--filter-expression "MQRankSum < -12.5" \
--filter-name "basic_snp_filter5" \
--filter-expression "ReadPosRankSum < -8.0" \
--filter-name "basic_snp_filter6" \
--filter-expression "SOR > 4.0"



gatk VariantFiltration \
-R genome.fa  \
-V raw_indels.vcf \
-O filtered_indels.vcf \
--filter-name "basic_snp_filter1" \
--filter-expression "QD < 2.0" \
--filter-name "basic_snp_filter2" \
--filter-expression "FS > 200.0" \
--filter-name "basic_snp_filter3" \
--filter-expression "ReadPosRankSum < -20.0" \
--filter-name "basic_snp_filter4" \
--filter-expression "SOR > 10.0" 
```
 
# Count variants 
```
grep 'PASS' filtered_indels.vcf | grep '1/1' -c
grep 'PASS' filtered_snps.vcf   | grep '1/1' -c
grep 'PASS' filtered_indels.vcf | grep '0/1' -c
grep 'PASS' filtered_snps.vcf   | grep '0/1' -c

grep '1/1' raw_indels.vcf -c
grep '1/1' raw_snps.vcf   -c
grep '0/1' raw_indels.vcf -c
grep '0/1' raw_snps.vcf   -c
```

