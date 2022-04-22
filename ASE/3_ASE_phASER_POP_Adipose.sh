#!/bin/bash


### This script generates haplotypic ASE data using phASER following procedure described in: https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
### Author: Shuli Liu shuliliu1991@yahoo.com 
### 4th. March, 2021

### Sould activate the conda environment
### conda activate ~/bin/conda_environment/env_python2.7/

###Dependencies: SciPy, NumPy, samtools, tabix, bedtools, Cython, pandas, IntervalTree,pysam

module load bcftools
module load samtools
module load bedtools
module load tabix

phaser="~/bin/phaser"
wk_dir="~/ASE"
gene_pos="~/ARS-UCD1.2/ARS1.2.genebody.bed"
vcf_dir="~/vcf_tissues/Adipose"
eQTL_dir="~/vcf_tissues/Adipose"
##################################
#### First remove all the _gene_ae.txt files to one folder.
cd ${wk_dir}
mkdir Gene_AE 
cp */*_phaser_gene_ae.txt ./Gene_AE

##############################
#### Construct the matrix of haplotypic expression: combine all individuals.

python ${phaser}/phaser_pop/phaser_expr_matrix.py --gene_ae_dir ./Gene_AE --features ${gene_pos} --o Haplotype_count_matrix

###################################
### Extract the combined vcf file of adipose samples

bcftools view -S Adipose.sample.txt -Oz ${vcf_dir}/Adipose.filtered.vcf.gz > Adipose.vcf.gz
tabix -fp vcf Adipose.vcf.gz

##################################
### generate the gene-snp pairs from eQTL results of Adipose

awk -F "_|\t" '{print $1,$2"_"$3"_"$4"_"$5, $2,$3,$4,$5}' OFS="\t" ${vcf_dir}/Adipose.nominals.sig.txt | sed '1i gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt' > Adipose.test_pairs.txt

python ${phaser}/phaser_pop/phaser_cis_var.py --bed Haplotype_count_matrix.gw_phased.bed.gz --vcf Adipose.vcf.gz --pair Adipose.test_pairs.txt --map Adipose_sample_map.txt --o Adipose_results.txt --ignore_v 1




