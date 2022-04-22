#!/bin/bash


##GWAS_TOOL Folder##
GWAS_TOOLS="~/bin/summary-gwas-imputation/src"

##GWAS File Folder##
Reference="~/run7_yak_indicus/parquet_run7"
OUTPUT="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS"

cd $OUTPUT/harmonized_gwas/
trait=`ls HM2_* | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed 's/HM2_//g; s/.csv.gz//g; s/.txt.gz//g' `

sh imputation.sh ${trait} 


