#!/bin/bash


########################################################################################################
# This script is used for testing the indicus specific cis-eQTLs for snp x breed interaction
# Shuli Liu shuliliu1991@cau.edu.cn shuliliu1991@yahoo.com 2021/07/04
#########################################################################################################
wk_dir="~/cis-eQTL_Muscle_test_SNP_breed/Muscle"
cd ${wk_dir}
#########Prepare the genotypes and phenotypes (adjust the covariances except the breed)
module load plink
module load vcftools
module load r
list="~/AF_tissue_specific_eQTLs/Muscle_indicus.specific.sig_pairs.txt" 
vcf="~/vcf_tissues/Muscle/Muscle.filtered.vcf.gz"

gene_list=`cat $list | cut -f 1| sort | uniq`

for gene in ${gene_list}
do
Rscript meta_analysis_cis-eQTL_prepare_phenotype.r Muscle ${gene}
grep ${gene} ${list} | cut -f 2 | uniq > ${gene}.snp_list.tmp
zcat ${vcf} | vcftools --vcf - --snps ${gene}.snp_list.tmp --recode --recode-INFO-all --out ${gene}.snp_only
plink --vcf ${gene}.snp_only.recode.vcf --pheno ${gene}.phenotype.txt.tmp --linear interaction --covar covariate.txt --chr-set 29 --allow-no-sex --parameters 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23  --out ${gene}
rm ${gene}.phenotype.txt.tmp ${gene}.snp_only.recode.vcf ${gene}.snp_list.tmp
done

grep "ADDxSpecies2" *.assoc.linear > ADDxSpecies_info.txt

Rscript meta_analysis_cis-eQTL_adjust_FDR.r