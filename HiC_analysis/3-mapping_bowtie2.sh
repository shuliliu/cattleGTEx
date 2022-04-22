#!/bin/bash

module load bwa
module load samtools
module load bowtie2
ref="~/Bovine_genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"

cd ~/Hi-C
mkdir mapping_bam
cd mapping_bam

bowtie2 -x  ${ref} --threads 16 -U ../trim_reads/*1_val_1.fq.gz --reorder --local| samtools view -Shb - > Lung_1_1.mapping.bam 



