#!/bin/bash


### This script generates haplotypic ASE data using phASER following procedure described in: https://stephanecastel.wordpress.com/2017/02/15/how-to-generate-ase-data-with-phaser/
### Author: Shuli Liu shuliliu1991@yahoo.com 
### 4th. March, 2021

### Sould activate the conda environment
### conda activate ~/conda_environment/env_python2.7/

###Dependencies: SciPy, NumPy, samtools, tabix, bedtools, Cython, pandas, IntervalTree,pysam

module load bcftools
module load samtools
module load bedtools
module load tabix

phaser="~/bin/phaser"
wk_dir="~/ASE/Adipose"
gene_pos="~/ASE/Adipose/Adipose_splicing.bed.txt"
vcf_dir="~/vcf_tissues/Adipose"
eQTL_dir="~/vcf_tissues/Adipose"
splicing="~/splicing_QTL"
##################################
#### First remove all the _gene_ae.txt files to one folder.
cd ${wk_dir}
mkdir Intron_AE 
cp */*_phaser_top_intron_ae.txt ./Intron_AE

##############################
#### Construct the matrix of haplotypic expression: combine all individuals.

python ${phaser}/phaser_pop/phaser_expr_matrix.py --gene_ae_dir ./Intron_AE --features ${gene_pos} --o Haplotype_count_matrix_intron --t 10

###################################
### Extract the combined vcf file of adipose samples

bcftools view -S Adipose.sample.txt -Oz ${vcf_dir}/Adipose.filtered.vcf.gz > Adipose.vcf.gz
tabix -fp vcf Adipose.vcf.gz

##################################
### generate the gene-snp pairs from eQTL results of Adipose

sed '1d' ${splicing}/Adipose_new1_perind2.permutation.storey.txt | awk -F "_|\t" '{print $1"_"$2"_"$3,$8"_"$9"_"$10"_"$11, $8,$9,$10,$11}' OFS="\t"  | sed '1i gene_id\tvar_id\tvar_contig\tvar_pos\tvar_ref\tvar_alt' > Adipose.test_pairs_intron.txt
sed -i 's/:/_/g' Adipose.test_pairs_intron.txt
zcat Haplotype_count_matrix_intron.gw_phased.bed.gz | sed 's/:/_/g' | uniq | bgzip -f > Haplotype_count_matrix_intron_formated.gw_phased.bed.gz
tabix -p bed -f Haplotype_count_matrix_intron_formated.gw_phased.bed.gz
python ${phaser}/phaser_pop/phaser_cis_var.py --bed Haplotype_count_matrix_intron_formated.gw_phased.bed.gz --vcf Adipose.vcf.gz --pair Adipose.test_pairs_intron.txt.sorted --map Adipose_sample_map.txt --o Adipose_results_intron.txt --ignore_v 1 --t 10




