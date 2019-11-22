---
title: 'How-to: allele specific expression' 
---

# Genotypes 
## Create a personal genome 
Using gtgtools to make a haploid version. Insert only the homozygous alternate SNPs/indels (ie markes as 1/1 in the VCF file). 


Other tools exist to do this too, such as WASP. This makes diploid genomes, and requires phased genomes.   



# Mapping 
See the RNA-seq alignment page.  

# Counting 
## IGVTools version 
```
igvtools=/sonas-hs/gillis/hpc/home/sballouz/IGVTools/igvtools.jar

java -jar $igvtools count -z 0 -w 1 --bases --strands read   \
   chrX.split.filtered.bam   \
   chrX.split.filtered.wig   \
   chrX.fa
```
## GATK
```
picard=/sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/picard.jar
GATK=/sonas-hs/gillis/hpc/home/sballouz/GenomeAnalysisTK/GenomeAnalysisTK.jar

java -jar $GATK \
   -T ASEReadCounter \
   -R genome.fa \
   -o chrX.counts.csv \
   -I chrX.split.filtered.bam \
   -sites chrX.split.filtered.vcf \
   -U ALLOW_N_CIGAR_READS
```
 
# Assessing 


