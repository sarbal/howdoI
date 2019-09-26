---
title: 'How-to: download data from SRA'
---
References: 
https://www.michaelgerth.net/news--blog/how-to-efficiently-bulk-download-ngs-data-from-sequence-read-databases


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
https://asperasoft.com/
https://downloads.asperasoft.com/downloads
https://downloads.asperasoft.com/connect2/ 
 
https://ftp-private.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/
https://www.ebi.ac.uk/ena/browse/read-download

https://gist.github.com/mfansler/71f09c8b6c9a95ec4e759a8ffc488be3

## Install aspera 
```
sh aspera-connect-[version].sh
```
## Download file 
```
ascp -QT -l 300m -P33001 -i <aspera connect installation directory>/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:<file or files to download> <download location>
```


# FTP from ENA 
```
wget “ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR101/006/SRR1016916/*”
```


