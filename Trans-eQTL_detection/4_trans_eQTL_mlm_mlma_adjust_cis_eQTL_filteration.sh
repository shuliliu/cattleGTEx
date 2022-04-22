#!/bin/bash


module load r
tissue=$1
wk_dir="Trans_eQTL_fit_GRM_mlma"

cd ${wk_dir}
cd ${tissue}
cd output_mlm_corrected
cat *.assoc.mlma.tmp | sed '/^Gene/d' | sed '1i Gene\tCHR\tSNP\tPOS\tA1\tA2\tAF1\tBETA\tSE\tP' > ../All_association.${tissue}.mlma_corrected.eQTL.txt
#rm *.assoc.mlma.tmp
cd ../
gene_num=`cat gene_list.txt | wc -l`
snp_num=`cat ${tissue}_GRM_input.bim | wc -l`
Rscript trans_eQTL_mlm_mlma_adjust_cis_eQTL_filteration.r ${tissue} ${gene_num} ${snp_num}
echo ${tissue}
