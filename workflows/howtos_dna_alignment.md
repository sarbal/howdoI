---
title: 'How-to: align DNA-seq data '
---
# Aligning with BWA

References: 
https://www.broadinstitute.org/files/shared/mpg/nextgen2010/nextgen_li.pdf
https://genomics.sschmeier.com/ngs-mapping/#mapping-statistics

## Get BWA 
http://sourceforge.net/projects/bio-bwa/files/
http://bio-bwa.sourceforge.net/bwa.shtml
```{}
bunzip2 bwa-0.X.X.tar.bz2 
tar xvf bwa-0.X.X.tar
cd bwa-0.X.X
make
export PATH=$PATH:/path/to/bwa-0.X.X
source ~/.bashrc
```
##   Download data 
###  Reference genome 
```{}
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.p12.genome.fa.gz

```

## Index genome
```{}
bwa index -a bwtsw GRCh38.genome.fa
```

##  Mapping 
###  Mapping Short reads
####  Paired-end
```{}
bwa aln -t 4 GRCh38bwaidx reads1.fastq.gz >  reads1.bwa
bwa aln -t 4 GRCh38bwaidx reads2.fastq.gz >  reads2.bwa
bwa sampe GRCh38bwaidx reads1.bwa reads2.bwa reads1.fastq.gz reads2.fastq.gz > aln.sam
```
####  Single-end
```{}
bwa aln -t 4 GRCh38bwaidx reads.fastq.gz >  reads.bwa
bwa samse GRCh38bwaidx reads.bwa reads.fastq.gz > aln.sam
```
### . Mapping long reads 
```{}
bwa bwasw GRCh38bwaidx long-reads.fastq.gz > aln-long.sam
```

### . Mapping with BWA-MEM
```{}
bwa index GRCh38.genome.fa
```
#### Paired-end
```{}
bwa mem GRCh38.genome.fa reads1.fastq.gz reads2.fastq.gz >  aln-pe.sam
```
####  Single-end
```{}
bwa mem GRCh38.genome.fa reads.fq.gz >  aln-se.sam
```


### Processing 
###  Fix mates and compress
```{}
samtools sort -n -O sam aln.sam | samtools fixmate -m -O bam - aln.fixmate.bam
```

###  Sorting
```{}
samtools sort -O bam -o aln.sorted.bam aln.fixmate.bam 
```

###  Remove duplicates
```{}
samtools markdup -r -S aln.sorted.bam aln.sorted.dedup.bam
```

