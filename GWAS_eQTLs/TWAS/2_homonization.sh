#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH -p short
#SBATCH -N 1 
#SBATCH -n 2
#SBATCH --time=48:00:00
#SBATCH --array=1-43 

##GWAS_TOOL Folder##
GWAS_TOOLS="~/bin/summary-gwas-imputation/src"

##GWAS File Folder##
DATA="~/Cattle_GWAS_ARS1.2"
Reference="~/run7_yak_indicus/parquet_run7"
OUTPUT="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS"
##Harmonization 1
cd ${OUTPUT}
gwas=`cat GWAS_sample_size.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 1`
n_sample=`cat GWAS_sample_size.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | cut -f 2`

python $GWAS_TOOLS/gwas_parsing.py \
-gwas_file $OUTPUT/harmonized_gwas/HM_${gwas} \
-snp_reference_metadata $Reference/variant_metadata.txt.gz METADATA \
-output_column_map variant_id variant_id \
-output_column_map non_effect_allele non_effect_allele \
-output_column_map effect_allele effect_allele \
-output_column_map effect_size effect_size \
-output_column_map pvalue pvalue \
-output_column_map chromosome chromosome \
--chromosome_format \
-output_column_map position position \
-output_column_map standard_error standard_error \
--insert_value sample_size ${n_sample} \
-output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size \
-output $OUTPUT/harmonized_gwas/HM2_${gwas}


