#!/bin/bash


#####################################################################################################
coloc_file="~/Colocalization/GWAS_hits"
LD_map="~/colocalization/LD_map"
GWAS="~/TWAS/2_RUN_SprediXcan/Prepare_GWAS/processed_summary_imputation"
eqtl="~/eQTL/cis-eQTL_Summary_Statistics"
wk_dir="~/colocalization/Example_all_plot_20211221"
TSS_list="~/ARS-UCD1.2/ARS1.2.TSS.bed"
######################################################################################################

cd ${wk_dir}

record=`head -n $SLURM_ARRAY_TASK_ID ${coloc_file}/enloc.rcp_0.5.sig_final_10-5_annot.txt | tail -n 1`
gene=`echo $record | awk -F ":|\t" '{print $1}' `
gene_name=`echo $record | awk '{print $11}'`
tissue=`echo $record | awk '{print $7}'`
trait=`echo $record | awk '{print $8}'`
target_snp=`echo $record | awk '{print $9}'`
chr=`echo ${target_snp} | cut -f 1 -d "_"`
pos=`grep $gene ${TSS_list} | cut -f 2`
start=$(($pos-1000000))
end=$(($pos+1000000))
target_pos=`echo ${target_snp} | cut -f 2 -d "_" `
rcp=`echo $record | awk '{print $6}'`



###Calculate LD matrix using plink.
module load bcftools
module load plink
module load r/4.1.2



zcat ${eqtl}/${tissue}.nominals.txt.gz | awk -v g=${gene} -F "_|\t| " '{if ($1==g) print $2, $3, -log($7)/log(10)}' OFS="\t" > ${trait}.${tissue}.${gene}.eqtl.txt
zat ${GWAS}/rsid_imputed_ARS_UCD1.2_formated.${trait}.*gz | awk -v a=chr${chr} -v s=${start} -v e=${end} '{s=s+0; e=e+0; if ($3==a &&  $4>=s && $4<=e ) print $3,$4,-log($11)/log(10)}' > ${trait}.${tissue}.${gene}.gwas.txt


###get the r2

plink --r2 --vcf ${LD_map}/Annot_Holstein.Chr${chr}-Run7-TAU-Beagle.vcf --ld-window 999999 --ld-snp ${target_snp} \
--ld-window-kb 1000 --ld-window-r2 0 --chr-set 29 --threads 4 --out ${trait}.${tissue}.${gene}.${target_snp}.r2.txt

awk '{print $1,$2,$3,$4,$5,$6,$7}' ${trait}.${tissue}.${gene}.${target_snp}.r2.txt.ld > ${trait}.${tissue}.${gene}.${target_snp}.r2.txt

awk 'FNR==NR{a[$5]=$7;next}{k=$2; if (k in a) print $0"\t"a[k]}' ${trait}.${tissue}.${gene}.${target_snp}.r2.txt  ${trait}.${tissue}.${gene}.eqtl.txt > ${trait}.${tissue}.${gene}.eqtl_r2.txt

awk 'FNR==NR{a[$5]=$7;next}{k=$2; if (k in a) print $0"\t"a[k]}' ${trait}.${tissue}.${gene}.${target_snp}.r2.txt  ${trait}.${tissue}.${gene}.gwas.txt > ${trait}.${tissue}.${gene}.gwas_r2.txt

##############################################
###Correlation between eqtl and gwas p-value.
############################################
Rscript locus_compare_plot.r ${gene} ${gene_name} ${tissue} ${trait} ${target_pos} ${target_snp} ${rcp}


