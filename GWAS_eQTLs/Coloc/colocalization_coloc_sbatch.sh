#!/bin/bash
#SBATCH --job-name=plot
#SBATCH -p short
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem 200G
#SBATCH -t 2-00:00:00

module load r

GWAS_dir="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation"
eQTL_dir="~/cis-eQTL"
wk_dir="coloc"
af_dir="~/vcf_tissues"

tissue=$1
gwas_file=$2
trait=`echo ${gwas_file} | cut -d "." -f 3`

eQTL_file=${eQTL_dir}/${tissue}.nominal.2rd.sig.maf_n_samples.txt

Rscript colocalization_coloc.r ${eQTL_file} ${wk_dir}/${trait}/${trait}.${tissue}.gwas.txt ${tissue} ${trait} ${wk_dir}/${trait}
