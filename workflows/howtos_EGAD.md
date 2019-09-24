---
title: "How-tos"
author: "Sara Ballouz"
date: "September 23, 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Alignment

## Aligning with STAR
1. Generating a genome index
a. Download files 
```{}

```

b. Run STAR to generate genome
```{}

```

2. Alignment 
a. One sample, paired end 
```{}
/sonas-hs/gillis/hpc/home/sballouz/STAR/STAR_2.7 \
 --genomeDir /sonas-hs/gillis/hpc/home/sballouz/sballouz/armadillo/ensembl/dasNov3.0.95/ \
 --readFilesIn /sonas-hs/gillis/hpc/home/sballouz/lyon/armadillo/test_fastq/12D101_062618_S3_R1_001.fastq.filtered_pairs_R1.fastq  /sonas-hs/gillis/hpc/home/sballouz/lyon/armadillo/test_fastq/12D101_062618_S3_R2_001.fastq.filtered_pairs_R2.fastq  \
 --readFilesCommand -  \
 --outSAMtype BAM Unsorted \
 --quantMode GeneCounts \
 --twopassMode Basic \
 --twopass1readsN -1
``` 

b. Multiple samples 

 
## Aligning with bowtie 


## Aligning with RSEM 
### Preparing reference:
```{}
rsem-prepare-reference --gtf ../GRCh38_Gencode22/gencode.v22.annotation.gtf ../GRCh38_Gencode22/GRCh38.p2.genome.fa.noPatches ../GRCh38_Gencode22_RSEM/RSEMgenome >& log &
```

### Quantification:
```{}
rsem-calculate-expression --bam --no-bam-output --seed 12345 -p 8 --paired-end --forward-prob 0 Aligned.toTranscriptome.out.bam ../GRCh38_Gencode22_RSEM/RSEMgenome Quant > & log.rsem
```

