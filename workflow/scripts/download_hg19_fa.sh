#!/bin/bash

mkdir ~/Downloads/hg19
cd hg19
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
mkdir uncompressed
for a in *.gz; do gunzip -c $a > uncompressed/`echo $a | sed s/.gz//`; done
cd uncompressed
cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrM.fa chrX.fa chrY.fa > hg19.fa

mkdir ~/Downloads/hg19_masked
cd hg19_masked
wget --timestamping 'http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFaMasked.tar.gz'
mkdir uncompressed_masked
for a in *.gz; do gunzip -c $a > uncompressed_masked/`echo $a | sed s/.gz//`; done
cd uncompressed_masked
cat chr1.fa.masked chr2.fa.masked chr3.fa.masked chr4.fa.masked chr5.fa.masked chr6.fa.masked chr7.fa.masked chr8.fa.masked chr9.fa.masked chr10.fa.masked chr11.fa.masked chr12.fa.masked chr13.fa.masked chr14.fa.masked chr15.fa.masked chr16.fa.masked chr17.fa.masked chr18.fa.masked chr19.fa.masked chr20.fa.masked chr21.fa.masked chr22.fa.masked chrM.fa.masked chrX.fa.masked chrY.fa.masked > hg19.masked.fa
