#!/bin/bash


module load bedtools
module load gsl


fastenloc="~/bin/fastenloc/src/fastenloc.static"
eQTL_annot="~/colocalization/Cattle_GTEx_eQTL"
gwas_dir="~/colocalization/GWAS_PIP"
out_dir="~/colocalization/output_enloc"

###We use gwas as a array factor. 
cd ${gwas_dir}

#gwas=`ls rsid_imputed_ARS_UCD1.2_formated.*.pip.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 | sed 's/.pip.gz//' `

gwas="rsid_imputed_ARS_UCD1.2_formated.Udder_cleft.G.EMAT" 

#######################################

total_variants=`zcat ${gwas}.zval.gz | wc -l`

cd ${eQTL_annot}
tissue=`ls *.vcf.gz | head -n $SLURM_ARRAY_TASK_ID | tail -n 1| cut -d "." -f 1 `
#for tissue in `ls *.vcf.gz | cut -d "." -f 1`
#
${fastenloc} -eqtl ${eQTL_annot}/${tissue}.fastenloc.eqtl.annotation.vcf.gz  -gwas ${gwas_dir}/${gwas}.pip.gz -total_variants ${total_variants} -t ${tissue} -thread 4 -prefix ${out_dir}/${tissue}.${gwas} -s 1
#done




