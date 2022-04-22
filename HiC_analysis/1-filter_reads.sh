#!/bin/bash

module load cutadapt

trim_galore="~/bin/TrimGalore-master/trim_galore"
cd ~/Hi-C/fastq

zcat *_1.fastq.gz | gzip -c > Lung_1.fq.gz
zcat *_2.fastq.gz | gzip -c > Lung_2.fq.gz

cd ../
mkdir trim_reads
cd trim_reads
${trim_galore} --paired --fastqc -o ./ ../fastq/Lung_1.fq.gz ../fastq/Lung_2.fq.gz

