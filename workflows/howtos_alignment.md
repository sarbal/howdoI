# How-to: align RNA-seq data 

# Aligning with STAR
## Get the latest version of STAR. 
All things STAR here: https://github.com/alexdobin/STAR. 

Precompiled binaries are easiest. 
```{}
wget https://github.com/alexdobin/STAR/blob/master/bin/Linux_x86_64_static/STAR
```

Or you can compile from source. 
```{}
# from latest release
wget https://github.com/alexdobin/STAR/archive/2.7.2b.tar.gz
tar -xzf 2.7.2b.tar.gz
cd STAR-2.7.2b

# or clone using git
git clone https://github.com/alexdobin/STAR.git

# compile (linux)
cd STAR-2.7.2b/source
make STAR
```


## Generating a genome index.
You will need a genome file (FASTA) and its corresponding annotation file (GTF). There are many locations and ways to download these files. GENCODE (https://www.gencodegenes.org/human/) has data for humans and mouse only, while ENSEMBL (http://useast.ensembl.org/index.html) and NCBI (https://www.ncbi.nlm.nih.gov/search/all/?term=human) will have other species. Note, these are also specific to the latest release and will change over time. 

- GRChXX: genome reference version
- pXX: patch version 
- GencodeXX: gene annotation version

### Download files. 
#### Genome (FASTA) file 
Pick one of these to download. 
```{}
mkdir GRCh38_Gencode31
cd GRCh38_Gencode31

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.p12.genome.fa.gz

wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz # unmasked
wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.fa.gz # masked with Ns
wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz # masked with lowercases ***

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/Assembled_chromosomes/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
```

#### GTF file with coordinates that correspond to genome files above

```{}
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz

wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GFF/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz
```

#### EntrezIDs

```{}
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
zgrep 9606 Homo_sapiens.gene_info.gz  | cut -f2,3 > 9606.ID2name
```

#### Barcodes from cellranger for 10x.
https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest. 
Need to download cellranger bundle to get to all the files. Note, need to go to website to get "permissions". 

```{}
wget http://cf.10xgenomics.com/releases/cell-exp/cellranger-3.1.0.tar.gz 
tar xvzf cellranger-3.1.0.tar.gz 
ls  cellranger-3.1.2/cellranger-cs/3.1.2/lib/python/cellranger/barcodes/
# 3M-february-2018.txt.gz  737K-april-2014_rc.txt  737K-august-2016.txt
```

Can also get genome + gtf files. 

```{}
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
```

#### Transcriptome (for kallisto)
```{}
wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
```

 
### Run STAR to generate genome. 
Run from the same directory or from the genome location. Specify output directories. Note, depending where you run this, you might need to modify RAM and number of threads. Overhang is default. Many other options can be tweaked and depend on the species and ultimate application. 
```{}
./STAR_2.7  \
--runThreadN 4       \
--runMode genomeGenerate  \
--genomeDir GRCh38_Gencode31   \
--genomeFastaFiles  GRCh38.p12.genome.fa  \
--sjdbGTFfile gencode.v31.annotation.gtf  \
--sjdbOverhang 100                 \
--limitGenomeGenerateRAM 40048000000
```

### Generate gene annotation file. 
This is helpful for downstream analyses. Note, ensemblIDs are not unique for XY paralogs and some tRNAs/rRNAs. This ignores them but maybe shouldn't. Needs some manual inspection when/if file is not exactly consistent (i.e., other species, new versions.)
```{r, eval=FALSE}
temp = read.table("gencode.v31.annotation.gtf.gz", header=F, sep="\t")
uid2symbol = read.table("9606.ID2name", header=F)
uid2symbol = unique(uid2symbol)

filt = temp[,3] == "gene"
mat = strsplit( as.character(temp[filt,9]), ";")

ids = sapply(1:length(mat), function(i) mat[[i]][1] )
genetype  = sapply(1:length(mat), function(i) mat[[i]][2] )
genestat  = sapply(1:length(mat), function(i) mat[[i]][3] )
genename  = sapply(1:length(mat), function(i) mat[[i]][3] )


ids = gsub( "gene_id", "", ids)
genetype = gsub( "gene_type", "", genetype)
genestat = gsub( "gene_status", "", genestat)
genename = gsub( "gene_name", "", genename)

ids = gsub( " ", "", ids)
genetype = gsub( "  ", "", genetype)
genestat = gsub( "  ", "", genestat)
genename = gsub( "  ", "", genename)


idssplit = matrix( unlist( strsplit(as.character(ids), "\\.") ), ncol=2, nrow=length(ids), byrow=T )
ids2 = gsub(" ", "", idssplit[,1] )

part1 = cbind( temp[filt,1:8], ids2, genetype, genename)

# EntrezIDs not provided in GTF file, need to add it 
m = match(uid2symbol[,2], part1[,12])
f.u = !is.na(m)
f.p = m[f.u]
part1 = cbind(part1,part1[,1])
part1[,13] = NA
part1[f.p,13] = uid2symbol[f.u,1]

# File format is not always consistent, need to be careful here, but mainly removing things/columns we don't need
attr = part1[,-2]
attr = attr[,-2]
attr = attr[,-4]

colnames(attr) = c("chr", "start","end", "strand","un", "ensemblID", "type", "stat", "name", "entrezID")
save(attr, file="gene_annotations_v31.Rdata")

```



## Alignment examples 
### Paired-end sample
#### Unsorted, BAM output
```{}
./STAR_2.7 \
 --genomeDir GRCh38_Gencode31  \
 --readFilesIn  sample_R1.fastq  sample_R2.fastq  \
 --readFilesCommand -  \
 --outSAMtype BAM Unsorted \
 --quantMode GeneCounts \
 --twopassMode Basic \
 --twopass1readsN -1
``` 

#### Multiple samples at once
```{}
./STAR_2.7 \
 --genomeDir GRCh38_Gencode31  \
 --readFilesIn  sample1_R1.fastq,sample2_R1.fastq,sample3_R1.fastq sample1_R2.fastq,sample2_R2.fastq,sample3_R2_001.fastq  \
 --readFilesCommand -  \
 --outSAMtype BAM Unsorted \
 --quantMode GeneCounts \
 --twopassMode Basic \
 --twopass1readsN -1
``` 

### Single-end samples
#### Mapping compressed files 
```{}
./STAR_2.7 \
 --genomeDir GRCh38_Gencode31   \
 --readFilesIn  single_end.fastq.gz    \
 --readFilesCommand zcat  \
 --quantMode GeneCounts \
 --twopassMode Basic \
 --twopass1readsN -1
``` 

#### Multiple samples
```{}
./STAR_2.7 \
 --genomeDir GRCh38_Gencode31   \
 --readFilesIn  single_end1.fastq.gz,single_end2.fastq.gz,single_end3.fastq.gz \
 --readFilesCommand zcat  \
 --quantMode GeneCounts \
 --twopassMode Basic \
 --twopass1readsN -1
``` 

### Single-cell samples from 10x 
Note, read 2 from output needs to come before read 1 because of the chemistry of 10x. If barcode length is not standard, it needs to be specified. 
```{}
./STAR_2.7  \
--runThreadN 1 \
--soloType Droplet \
--genomeDir GRCh38_Gencode31 \
--outSAMtype BAM Unsorted \
--quantMode GeneCounts \
--soloCBwhitelist 3M-february-2018.txt \
--readFilesIn sc1R2.fastq,sc2_R2.fastq,sc3_R2.fastq,sc4_R2.fastq sc1_R1.fastq,sc2_R1.fastq,sc3_R1.fastq,sc4_R1.fastq \
--soloBarcodeReadLength 28 
```



### Parsing STAR output in R
#### Example from multiple bulk runs
```{r, eval=FALSE}
files = as.character(unlist(read.table("runs") ))
genecounts = "ReadsPerGene.out.tab"
splicejunctions = "SJ.out.tab"
logname = "Log.final.out"
load("gene_annotations.Rdata")
dir = "outs/"

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

#### Example from single-cell data 
```{}
# Filter barcodes 
ls Solo.out 
barcodes.tsv  Gene.stats  genes.tsv  matrix.mtx
```

Adapted from: https://davetang.org/muse/2018/08/09/getting-started-with-cell-ranger/ 
```{r, eval=FALSE}
library(DropletUtils)
library(dplyr)
library(ggplot2)

# Load STAR solo output (reads in matrix.mtx file)
data <- read10xCounts(data.dir = "Solo.out/")

# Rank barcodes to find knee and inflection points -> where UMI counts shift suggesting empty drops/ambient RNA
br.out <- barcodeRanks(counts(sce))
br.out.df <- as.data.frame(br.out)
br.out.df$barcode <- colData(sce)$Barcode
br.out.df %>% filter(rank <= 10) %>% arrange(rank)
x_knee <- br.out.df %>% filter(total > br.out$knee) %>% arrange(total) %>% select(rank) %>% head(1)
x_inflection <- br.out.df %>% filter(total > br.out$inflection) %>% arrange(total) %>% select(rank) %>% head(1)
padding <- length(br.out$rank) / 10
  

# Bend the knee!  
ggplot(br.out.df, aes(x = rank, y = total)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        title = element_text(size = 16)) +
  geom_hline(yintercept = br.out$knee, linetype = 2, colour = "dodgerblue") +
  geom_hline(yintercept = br.out$inflection, linetype = 2, colour = "forestgreen") +
  labs(x = "Rank", y = "Total", title = "Barcode Rank vs Total UMI") +
  annotate("text", label = paste0("Knee (", x_knee, ")"), x = x_knee$rank + padding, y = br.out$knee, size = 5) +
  annotate("text", label = paste0("Inflection (", x_inflection, ")"), x = x_inflection$rank + padding, y = br.out$inflection, size = 5)


# Find the empty drops 
e.out <- emptyDrops(counts(sce))
 
# use FDR threshold of 0.01
is.cell <- e.out$FDR <= 0.01
 
# replace all NA's with FALSE
is.cell.no.na <- is.cell
is.cell.no.na[is.na(is.cell)] <- FALSE
sum(is.cell.no.na)
 

# Plot 
ggplot(br.out.df[is.cell.no.na,], aes(x = rank, y = total)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        title = element_text(size = 16)) +
  labs(x = "Rank", y = "Total", title = "Cell barcodes that differ from ambient RNA")

 # 
 sce <- sce[,is.cell.no.na]
 save(sce, file="sce.Rdata")

```


# Aligning with cellranger 
Since all the necessary referenes come with cellranger (at least human), this is straightforward. For other species, one can repeat the STAR genome index generation. Here is an example from a PBMC human run. 
The fastq files from the sequencing machine should have read 1, read 2 and an index file.  
```{}
ls fastq_path
Gillis_05_XCGD10xGEX_PBMC_S1_L001_I1_001.fastq.gz
Gillis_05_XCGD10xGEX_PBMC_S1_L001_R1_001.fastq.gz  
Gillis_05_XCGD10xGEX_PBMC_S1_L001_R2_001.fastq.gz  
Gillis_05_XCGD10xGEX_PBMC_S1_L002_I1_001.fastq.gz  
Gillis_05_XCGD10xGEX_PBMC_S1_L002_R1_001.fastq.gz 
Gillis_05_XCGD10xGEX_PBMC_S1_L002_R2_001.fastq.gz  
Gillis_05_XCGD10xGEX_PBMC_S1_L003_I1_001.fastq.gz  
Gillis_05_XCGD10xGEX_PBMC_S1_L003_R1_001.fastq.gz  
Gillis_05_XCGD10xGEX_PBMC_S1_L003_R2_001.fastq.gz
Gillis_05_XCGD10xGEX_PBMC_S1_L004_I1_001.fastq.gz
Gillis_05_XCGD10xGEX_PBMC_S1_L004_R1_001.fastq.gz
Gillis_05_XCGD10xGEX_PBMC_S1_L004_R2_001.fastq.gz
```
Cellranger maps and quantifies. One giant bam file is typically produced. This can be parsed for further analyses, and also is a way to capture the experiment as pseudobulk.
```{}
cellranger count --id=Gillis_PBMC \
        --jobmode=sge \
        --transcriptome=refdata-cellranger-GRCh38-3.0.0 \
        --fastqs=fastq_path \
        --sample=Gillis_PBMC \
        --project=Gillis
```


# Aligning with RSEM 

STAR outputs count summaries only. If you wish for quantified data (primarily bulk), either summarized to the transcripts or to genes, you can run RSEM. It can take in the bam files from STAR (or remap), and outputs TPMs. 
## Get RSEM 
https://deweylab.github.io/RSEM/ 
https://github.com/bli25broad/RSEM_tutorial
```{}
git clone git@github.com:bli25ucb/RSEM_tutorial.git 

cd software
unzip bowtie2-2.2.6-source.zip

cd bowtie2-2.2.6
make -j 8

cd ..
tar -xzf RSEM-1.2.25.tar.gz

cd RSEM-1.2.25
make -j 8
make ebseq
cd ..
```


## Preparing the reference

```{}
mkdir GRCh38_Gencode31_RSEM 

rsem-prepare-reference \ 
   --gtf  GRCh38_Gencode31/gencode.v31.annotation.gtf \ 
   GRCh38_Gencode31/GRCh38.p12.genome.fa  GRCh38_Gencode31_RSEM/RSEMgenome >& log &
```

## Alignment (with bowtie)  
```{}
rsem-calculate-expression -p 8 --paired-end \
					--bowtie2 --bowtie2-path  bowtie2 \
					--estimate-rspd \
					--append-names \
					--output-genome-bam \
					sample_reads1.fastq sample_reads2.fastq \
					GRCh38_Gencode31_RSEM/RSEMgenome experimentID
```
 

## Quantification (from mapped)
### Paired-end 
```{}
rsem-calculate-expression \ 
   --bam --no-bam-output \ 
    -p 8 --paired-end --forward-prob 0 \ 
    output/Aligned.toTranscriptome.out.bam \
    GRCh38_Gencode31_RSEM/RSEMgenome experimentID 
```
 
 
# Aligning with kallisto 
Pseudoalignments. 

## Get kallisto
https://pachterlab.github.io/kallisto/download.html

## Preparing the reference
Kallisto uses the transcriptome rather than the genome. Pre-indexed versions here:
https://github.com/pachterlab/kallisto-transcriptome-indices
Alternatively:

```{}
kallisto index -i GRCh38.cdna.idx GRCh38.cdna.all.fa.gz
```

##  Quantification

###  Paired-end
```{}
kallisto quant -i GRCh38.p12.genome.idx -o output -b 100 sample_R1.fastq sample_R2.fastq
```

###  Single-end

```{}
kallisto quant -i GRCh38.p12.genome.idx -o output -b 100 --single -l 180 -s 20 sample.fastq
```

###  Single-cell version 
https://bustools.github.io/BUS_notebooks_R/10xv2.html 


# Aligning with Salmon

## Get salmon
Download from: https://github.com/COMBINE-lab/salmon
Need some of the tools from: https://github.com/COMBINE-lab/SalmonTools
Pre-computed decoy files: https://drive.google.com/drive/folders/14VqSdZAKH82QwDWhMXNLFqMoskoqv3fS
e.g., human: 

```{}
https://drive.google.com/open?id=1T0_sQxCXbnLfT3oo764eKG2l743swn6e
```

##  Preparing the reference 
Like kallisto, salmon uses the transcriptome. 

```{}
GRCh38.cdna.all.fa.gz

salmon index -t GRCh38.cdna.all.fa -i transcripts_index -decoys decoys.txt -k 31 
```


## Alignment and quantification

###  Paired-end

####  Automatic library detection

```{}
salmon quant -i transcripts_index -l A -1 sample_R1.fastq -2 sample_R2.fastq --validateMappings -o transcripts_quant
```

####  Specified library type 

```{}
salmon quant -i transcripts_index -l IS -1 sample_R1.fastq -2 sample_R2.fastq --validateMappings -o transcripts_quant
```

####  Multiple files 

```{}
salmon quant -i transcripts_index -l A -1 sample1_R1.fastq sample2_R1.fastq -2 sample1_R2.fastq sample2_R2.fastq --validateMappings -o transcripts_quant
```

###  Single-end

```{}
salmon quant -i transcripts_index -l A -r sample.fastq --validateMappings -o transcripts_quant
```


##  Quantification

```{}
salmon quant -t transcripts.fa -l A -a aln.bam -o salmon_quant
```


###  Single-cell version 
https://combine-lab.github.io/alevin-tutorial/#blog


