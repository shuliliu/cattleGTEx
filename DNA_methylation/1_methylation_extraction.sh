#!/bin/bash

module load fastqc
module load bowtie2
module load cutadapt
module load samtools
module load sratoolkit


file=`ls SRR552*| head -n $SLURM_ARRAY_TASK_ID|tail -n 1`


name=`echo ${file}|cut -d "." -f 1`
trim_galore="~/bin/TrimGalore-master/trim_galore"
bismark="~/bin/Bismark_v0.19.0"
genome_folder="~/dor_new_assembly"

#run trim_galore

mkdir ${name}
cd ${name}
fastq-dump --split-3 --gzip ~/${file}
mkdir trim_reads

${trim_galore} --paired --fastqc --max_n 15 -o ./trim_reads ${name}_1.fastq.gz ${name}_2.fastq.gz


#run bismark
mkdir bamfile
cd ./trim_reads

if [ -e ${name}_1_val_1.fq.gz ]
then
${bismark}/bismark --multicore 2 --bowtie2 --gzip -p 4 -N 0 -o ../bamfile ${genome_folder} -1 ${name}_1_val_1.fq.gz -2 ${name}_2_val_2.fq.gz
else
exit
fi


#run deduplicate_bismark
cd ../bamfile

if [ -e ${name}_1_val_1_bismark_bt2_pe.bam ]
then
${bismark}/deduplicate_bismark -p --bam ${name}_1_val_1_bismark_bt2_pe.bam
else
exit
fi

mkdir ../methylation
${bismark}/bismark_methylation_extractor -p --gzip --ignore_r2 6 --multicore 8 --bedgraph -o ../methylation --cytosine_report --genome_folder ${genome_folder} *.deduplicated.bam



             
