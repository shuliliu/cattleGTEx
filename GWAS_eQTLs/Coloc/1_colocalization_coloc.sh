#!/bin/bash


module load r

GWAS_dir=~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation
eQTL_dir=~/cis-eQTL
wk_dir=coloc


cd ${GWAS_dir}

for gwas_file in `ls rsid_imputed_ARS_UCD1.2_formated.*`
do

trait=`echo ${gwas_file} | cut -d . -f 3`
mkdir ${wk_dir}/${trait}

cd ${wk_dir}
for tissue in Adipose Blood Embryo Hypothalamus Ileum Intramuscular_fat Jejunum Leukocyte Liver Lung Lymph_node Macrophage Mammary Milk_cell Monocytes Muscle Ovary Oviduct Pituitary Rumen Salivary_gland Skin_fibroblast Testis Uterus
do

sbatch colocalization_coloc_sbatch.sh ${tissue} ${gwas_file}
done

done
