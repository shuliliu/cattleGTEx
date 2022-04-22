#!/bin/bash


module load samtools
module load r
module load bowtie2
cd ~/Hi-C

~/bin/Software/HiC-Pro_2.11.4/bin/HiC-Pro -i ~/Hi-C/fastq -o ~/Hi-C/OUTPUT_HiC-Pro -c ~/Hi-C/config-hicpro.txt -p