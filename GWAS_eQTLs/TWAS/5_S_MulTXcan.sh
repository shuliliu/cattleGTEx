#!/bin/bash




##############################################################################
##Run MultiXcan to aggregated evidence across all tissues with S-multiXcan. 
#############################################################################
models="~/TWAS/PredictDB_Pipeline_GTEx_v7/model_training"
result_dir="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/TWAS_results"
GWAS_dir="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation"
output="."

source ~/.bashrc
source ~/.bash_profile

conda activate ~/bin/conda_environment/env_python3.7/


cd ${GWAS_dir}
gwas=`ls rsid_imputed_ARS_UCD1.2_formated*.txt.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed 's/rsid_//; s/.gz//' `
trait=`echo ${gwas} | awk -F "." '{print $3}'`


metaXcan="~/TWAS/MetaXcan/software"
python ${metaXcan}/SMulTiXcan.py \
--models_folder ${models}/dbs \
--models_name_pattern "gtex_v1_(.*)_imputed_Bos_signif.db" \
--models_name_filter "gtex_v1_(.*)_imputed_Bos_signif.db" \
--snp_covariance ~/TWAS/PredictDB_Pipeline_GTEx_v7/model_training/smultixcan_covariances.txt.gz \
--metaxcan_folder  ${result_dir}/${trait} \
--metaxcan_filter "twas_imputed_ARS_UCD1.2_formated.${trait}.G.(.*).EMAT.txt" \
--metaxcan_file_name_parse_pattern "twas_imputed_ARS_UCD1.2_formated.(.*).G.(.*).EMAT.txt" \
--gwas_file rsid_${gwas}.gz \
--snp_column variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--pvalue_column pvalue \
--keep_non_rsid \
--model_db_snp_key varID \
--cutoff_condition_number 30 \
--verbosity 7 \
--output ${output}/${gwas}_smultiXcan.txt \
--throw




