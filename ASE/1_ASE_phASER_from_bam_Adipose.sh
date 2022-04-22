#!/bin/bash


###

### This script generates haplotypic ASE data using phASER following procedure described in: https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
### Author: Shuli Liu shuliliu1991@yahoo.com 
### 3th. March, 2021

### Sould activate the conda environment
### conda activate ~/bin/conda_environment/env_python2.7/

###Dependencies: SciPy, NumPy, samtools, tabix, bedtools, Cython, pandas, IntervalTree,pysam

module load bcftools
module load samtools
module load bedtools
module load tabix

###FOR Adipose as an example
wk_dir="~/ASE/Adipose"
vcf_dir="~/vcf_tissues/Adipose"
bam_dir="~/CAREFUL_bam_files"
phaser="~/bin/phaser"
mappability="~/Cow_alignability"
gene_pos="~/ARS-UCD1.2/ARS1.2.genebody.bed"
cd ${wk_dir}

sample=`cat Adipose.deduplicated.sample.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`

reads=`cat Adipose.deduplicated.sample.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`

mkdir ${sample}
cd ${sample}

if [ ${reads} == "SINGLE" ]
then
    yes=0
else
    yes=1
fi
#####################################
#extract vcf file from combined vcf file
bcftools view -s ${sample} -Oz ${vcf_dir}/Adipose.filtered.vcf.gz > ${sample}.Adipose.vcf.gz

tabix -fp vcf ${sample}.Adipose.vcf.gz
###################################
#index bam file
samtools index ${bam_dir}/*/${sample}-STARAligned.sortedByCoord.out.bam

###################################
### Prepare the bed file with mappability < 0.5

awk '{if ($4<0.5) print}' ${mappability}/bosTau9_mappability_75mer_2mismatch.bed  > bosTau9_mappability_75mer_2mismatch_0.5.bed

python ${phaser}/phaser/phaser.py --vcf ${sample}.Adipose.vcf.gz --bam ${bam_dir}/*/${sample}-STARAligned.sortedByCoord.out.bam --paired_end ${yes} --mapq 255 --baseq 10 --sample ${sample} --haplo_count_blacklist bosTau9_mappability_75mer_2mismatch_0.5.bed  \
--threads 4 --o ${sample}_phaser \
--gw_phase_vcf 1


####################################
### Generate haplotype expression quantifications

python ${phaser}/phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts ${sample}_phaser.haplotypic_counts.txt --features ${gene_pos} --o ${sample}_phaser_gene_ae.txt