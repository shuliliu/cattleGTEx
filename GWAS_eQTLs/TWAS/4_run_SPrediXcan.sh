#!/bin/bash



SPrediXcan="~/TWAS/MetaXcan/software/SPrediXcan.py"

GWAS_dir="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation"
model_dir="~/TWAS/PredictDB_Pipeline_GTEx_v7/model_training"
result_dir="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/TWAS_results"

cd ${GWAS_dir}
gwas=`ls imputed_* | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
cd ${result_dir}
cd ${GWAS_dir}

for tissue in `cat ${model_dir}/scripts/gtex_tissues.txt`
do
python ${SPrediXcan} \
--model_db_path ${model_dir}/dbs/gtex_v1_${tissue}_imputed_Bos_signif.db \
--covariance ${model_dir}/covariances/gtex_v1_${tissue}_imputed_Bos_covariances.txt.gz \
--gwas_file ${gwas} \
--snp_column variant_id \
--effect_allele_column effect_allele \
--non_effect_allele_column non_effect_allele \
--zscore_column zscore \
--pvalue_column pvalue \
--keep_non_rsid \
--overwrite \
--output_file ${result_dir}/${tissue}/twas_${gwas}
done

