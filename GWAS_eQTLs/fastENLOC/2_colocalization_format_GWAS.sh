#!/bin/bash


module load bedtools
#### Require imputed GWAS summary statistics: Z-score
#### for Holstein GWAS statistics: Haplotype block map: use Holstein, constructed from 866 animals
#### columns are: ID, Block number, Z-score


LD_map="~/LD_map/HB_holstein.blocks.total.txt"
imputed_GWAS="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation"
wk_dir="~/colocalization/GWAS_PIP"
torus="~/bin/torus/src/torus.static"

###format the LD_map into chr BP1 BP2 BLOCK_ID
cd ${wk_dir}
#sed '1d' ${LD_map} | awk '{print $1,$2,$3, "LOC"NR}' OFS="\t" > HB_holstein.bed  ###LD map is only from holstein.
cd ${imputed_GWAS}
gwas=`ls rsid_imputed_ARS_UCD1.2_formated.*.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
name=`echo ${gwas} | sed 's/.txt.gz//; s/.csv.gz//' `
zcat ${gwas} | sed '1d' | awk '{print $3,$4,$4,$10,$1}' OFS="\t" | sed 's/chr//g' | bedtools intersect -a - -b ${wk_dir}/HB_holstein.bed -wa -wb | awk '{print $5,$9,$4}' OFS="\t" | gzip -c > ${wk_dir}/${name}.zval.gz
${torus} -d ${wk_dir}/${name}.zval.gz --load_zval -dump_pip ${wk_dir}/${name}.pip
gzip -f ${wk_dir}/${name}.pip





