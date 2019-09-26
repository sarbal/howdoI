---
title: 'How-to: download data from SRA'
---

# Using SRA toolkit
https://ncbi.github.io/sra-tools/install_config.html

## Install sratoolkit 
```
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz"
tar -xzf sratoolkit.current-centos_linux64.tar.gz
```
## Configure
```
./vdb-config -i
```
Follow the instructions to change the download directory. The tools will download the data into /ncbi/public/sra/, wherever you have it set. When downloading protected data such as GTEx, you will need a project key (prj_XXXX.ngc) and to configure the toolkit to use it. 

## Download a sample as an SRA file 
```
prefetch SRR3485764
fastq-dump –split-files /ncbi/public/sra/SRR3485764.sra 
```

## Download a sample as a fastq file 
```
fastq-dump –split-files SRR3485764 
```
Note, this still downloads the .sra file first. 

# Using aspera connect 
## Install aspera 


# Direct downloads 
## FTP from SRA 
## FTP from EBI 
