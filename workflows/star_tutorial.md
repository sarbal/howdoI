---
title: 'The STARs align'
---
 
# Setting up 
Hopefully, everyone can log into their (or shared) rugen node. Make a directory for this workshop. 
```{}
mkdir STAR_tutorial
cd STAR_tutorial
```

# Install STAR 
Precompiled binaries are easiest. 
```{}
wget https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64_static/STAR

wget https://github.com/alexdobin/STAR/blob/master/bin/MacOSX_x86_64/STAR
```

# Get data 
You will need a genome file (FASTA) and its corresponding annotation file (GTF). There are many locations and ways to download these files. GENCODE (https://www.gencodegenes.org/human/) has data for humans and mouse only, while ENSEMBL (http://useast.ensembl.org/index.html) and NCBI (https://www.ncbi.nlm.nih.gov/search/all/?term=human) will have other species. Note, these are also specific to the latest release and will change over time. 

- GRChXX: genome reference version
- pXX: patch version 
- GencodeXX: gene annotation version

## Genome (FASTA) and annotation (GTF) file 
```{}
mkdir GRCh38_Gencode31
cd GRCh38_Gencode31
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.p12.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
gunzip * 
```

## Reads/experiment
I have two experiments, but they are a little large and we wont have time to wait for the alignment results during this workshop. 
But, the inputs (and outputs) can be found in this folder on tyrone:  
```{}
ls /data/sballouz
ase_run  bulk_run  single_cell_run
```
Feel free to copy or link. 

# Prepare genome 
Run from the same directory or from the genome location. 
Specify output directories. 
Note, depending where you run this, you might need to modify RAM and number of threads. 
Overhang is default. Many other options can be tweaked and depend on the species and application. 
```{}
STAR  \
--runThreadN 10       \
--runMode genomeGenerate  \
--genomeDir GRCh38_Gencode31   \
--genomeFastaFiles  GRCh38.p12.genome.fa  \
--sjdbGTFfile gencode.v31.annotation.gtf  \
--sjdbOverhang 100                 \
--limitGenomeGenerateRAM 40048000000
```

# Alignment 
## Bulk RNA-seq experiment
```{}
STAR \
--genomeDir /data/genomes/GRCh38_Gencode31/ \
--readFilesIn /data/sballouz/bulk_run/xcgd_carrier_1_pbmc_S4_R1_001.fastq  /data/sballouz/bulk_run/xcgd_carrier_1_pbmc_S4_R2_001.fastq \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--twopassMode Basic \
--twopass1readsN -1 \
--runThreadN 10
```

## Single-cell RNA-seq experiment 
```{}
STAR \
--runThreadN 20 \
--soloType Droplet \
--genomeDir /data/genomes/GRCh38_Gencode31/ \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts \
--soloCBwhitelist /data/genomes/cellranger_barcodes/3M-february-2018.txt \
--readFilesIn SCGC-GILL-JG-03_S1_L001_R2_001.fastq.gz,SCGC-GILL-JG-03_S1_L002_R2_001.fastq.gz,SCGC-GILL-JG-03_S1_L003_R2_001.fastq.gz,SCGC-GILL-JG-03_S1_L004_R2_001.fastq.gz SCGC-GILL-JG-03_S1_L001_R1_001.fastq.gz,SCGC-GILL-JG-03_S1_L002_R1_001.fastq.gz,SCGC-GILL-JG-03_S1_L003_R1_001.fastq.gz,SCGC-GILL-JG-03_S1_L004_R1_001.fastq.gz \
--readFilesCommand zcat \
--soloCellFilter  CellRanger2.2 8000 0.99 10 \
--soloFeatures Gene Velocyto \
--soloBarcodeReadLength 28 
```

# Inspecting output 
Start R, and then read in the bulk sample output.  
```{r}
files = "STAR_out"
genecounts = "ReadsPerGene.out.tab"
splicejunctions = "SJ.out.tab"
logname = "Log.final.out"
# load("gene_annotations.Rdata")
dir = "."

Ns = list()
i = 1

for( n in files ){
  N = list()
  filedir = paste(dir, n, sep="/")
  countfile = paste(filedir, genecounts, sep=".")
  logfile = paste(filedir, logname, sep=".")
  if( file.exists(countfile) ) {
    print(countfile)
    
    counts =  read.table(countfile)
    log1 =read.table(logfile, sep="\t", nrows=6)
    log2 =read.table(logfile, sep="\t", skip=8, nrows=14)
    log3 =read.table(logfile, sep="\t", skip=23, nrows=4)
    log4 =read.table(logfile, sep="\t", skip=28, nrows=3)

    N$mapinfo      = rbind(log1,log2,log3,log4)
    N$unmapped     = counts[1,]
    N$multimapping = counts[2,]
    N$noFeature    = counts[3,]
    N$ambiguous    = counts[4,]
    N$length       = dim(counts)[1]-4
    N$genes        = counts[(1:N$length)+4,1]
    N$counts1      = counts[(1:N$length)+4,2]
    N$counts2      = counts[(1:N$length)+4,3]
    N$counts3      = counts[(1:N$length)+4,4]
 } else {
   # Stranded or unstranded? Spot check this before running the code. Should use the last column if stranded. Otherwise first/max. 
   N$counts3 = rep(0, length(attr$ensemblID ) )
 }
 if( i > 1  ){
   counts_exp = cbind(counts_exp, N$counts3)
 } else {
   counts_exp = N$counts3
 }
 Ns[[i]] = N
 print(i)
 i = i + 1
}

rownames(counts_exp) = attr$ensemblID
colnames(counts_exp) = files
save(Ns, counts_exp, file=paste(dir, "counts.Rdata", sep=""))
```



