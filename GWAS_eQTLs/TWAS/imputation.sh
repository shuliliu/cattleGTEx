#!/bin/bash

date

##GWAS_TOOL Folder##
GWAS_TOOLS="~/bin/summary-gwas-imputation/src"

##GWAS File Folder##
Reference="~/RNA-seq/run7_yak_indicus/parquet_run7"
OUTPUT="~/RNA-seq/TWAS/2_RUN_SprediXcan/Prepare_GWAS"
LD="~/Bovine_genome"


trait=$1
for chr in $(seq 1 29)
do
for i in $(seq 0 9)
do
python $GWAS_TOOLS/gwas_summary_imputation.py \
-gwas_file $OUTPUT/harmonized_gwas/HM2_${trait}.* \
-by_region_file ${LD}/genome.1mb.win.bed.gz \
-parquet_genotype ${Reference}/ARS_UCD1.2__no_MAF_monoallelic_variants.chr${chr}.variants.parquet \
-parquet_genotype_metadata ${Reference}/ARS_UCD1.2__no_MAF_monoallelic_variants.variants_metadata.parquet \
-window 100000 \
-parsimony 7 \
-chromosome ${chr} \
-regularization 0.1 \
-sub_batches 10 \
-sub_batch ${i} \
--standardise_dosages \
-output $OUTPUT/summary_imputation/${trait}.chr${chr}_${i}.txt.gz
done
done

python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
-gwas_file $OUTPUT/harmonized_gwas/HM2_${trait}.*.gz \
-folder $OUTPUT/summary_imputation \
-pattern ${trait}.* \
-parsimony 7 \
-output $OUTPUT/processed_summary_imputation/imputed_${trait}.txt.gz
date



