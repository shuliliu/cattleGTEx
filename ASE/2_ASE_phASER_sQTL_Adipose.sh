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
splicing="~/cis-sQTL/Adipose"
cd ${wk_dir}

sample=`cat Adipose.deduplicated.sample.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`

reads=`cat Adipose.deduplicated.sample.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`

cd ${sample}


################################################
##Construct the bed file for intron regions.
#### cd ${wk_dir}
 awk -F ":|\t" '{print $1,$2,$3, $1":"$2":"$3":"$4}' OFS="\t" ${splicing}/Adipose.nominal.chrAll.sig.txt > Adipose_splicing.bed


####################################
### Generate haplotype expression quantifications

python ${phaser}/phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts ${sample}_phaser.haplotypic_counts.txt --features Adipose_splicing.bed --o ${sample}_phaser_top_intron_ae.txt
rm *.tmp
