---
title: 'How-to: call variants from RNA-seq data'
---
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11164


# Install GATK4 
Download the precompiled jar binaries. Set paths. 

# Adding read groups to bam file 
```
java -jar picard AddOrReplaceReadGroups \
    I=sample.filt.bam \
    O=sample.rg.bam \
    SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
```

# Marking duplicates
```
java -jar picard MarkDuplicates \
    I=sample.rg.bam  \
    O=sample.dedupped.bam \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
```

# Splitting and trimming
```
gatk SplitNCigarReads \
   -R genome.fa \
   -I sample.dedupped.bam \
   -O sample.split.filtered.bam  
```

# Haplotype calling  
```
gatk HaplotypeCaller \
   -R genome.fa \
   -L 1 2 3 4 5 6 X \
   -I  sample.split.filtered.bam \
   -O  sample.filtered.vcf \
   --dont-use-soft-clipped-bases \
   --standard-min-confidence-threshold-for-calling 20.0 
```


