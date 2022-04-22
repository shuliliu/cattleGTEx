#!/bin/bash



##############################Trans-eQTL fit GRM mlm######################################
##Author Shuli Liu shuliliu1991@yahoo.com
##Date: 2021.03.24
##########################################################################################

##load required software
module load plink
module load r
## The required directories and files
tissue=$1
used_genes="~/ARS-UCD1.2/ARS1.2.genebody.protein_coding-lncRNA.bed" ##23341 genes
wk_dir="Trans_eQTL_fit_GRM_mlma"

cd ${wk_dir}
cd ${tissue}

gene_file=`ls gene_list*.tmp | head -n $SLURM_ARRAY_TASK_ID | tail -n 1 `

gene_list=(`cat ${gene_file}`)
for i in "${!gene_list[@]}"
do
num_batch=$((SLURM_ARRAY_TASK_ID-1))
pheno_id=$((i+300*num_batch+1))
if [ -e fixed.txt ]; then
gcta64 --mlma --bfile ${tissue}_GRM_input --grm ${tissue}_geno_grm --autosome-num 29 \
--pheno phenotype.txt --mpheno ${pheno_id}  --qcovar cvrt.txt --covar fixed.txt  --reml-maxit 2000 --reml-bendV \
--thread-num 10 --out output_mlm/${gene_list[i]}.assoc
if [ $? -gt 0 ]; then
echo "Failed attempt of gcta gwas of ${gene_list[i]} with full ov and reml alg 0"
gcta64 --mlma --bfile ${tissue}_GRM_input --grm ${tissue}_geno_grm --autosome-num 29 --reml-alg 1 \
--pheno phenotype.txt --mpheno ${pheno_id}  --qcovar cvrt.txt --covar fixed.txt  --reml-maxit 2000 --reml-bendV  \
--thread-num 10 --out output_mlm/${gene_list[i]}.assoc
if [ $? -gt 0 ]; then
echo "Failed attempt of gcta gwas of ${gene_list[i]} with reml alg 1"
gcta64 --mlma --bfile ${tissue}_GRM_input --grm ${tissue}_geno_grm --autosome-num 29 --reml-alg 2 \
--pheno phenotype.txt --mpheno ${pheno_id}  --qcovar cvrt.txt --covar fixed.txt  --reml-maxit 2000 --reml-bendV  \
--thread-num 10 --out output_mlm/${gene_list[i]}.assoc
if [ $? -gt 0 ]; then
echo "Failed attempt of gcta gwas of ${gene_list[i]} with reml alg 2"
fi
fi
fi

else
gcta64 --mlma --bfile ${tissue}_GRM_input --grm ${tissue}_geno_grm --autosome-num 29 \
--pheno phenotype.txt --mpheno ${pheno_id} --qcovar cvrt.txt  --reml-maxit 2000 --reml-bendV \
--thread-num 10 --out output_mlm/${gene_list[i]}.assoc
if [ $? -gt 0 ]; then
echo "Failed attempt of gcta gwas of ${gene_list[i]} with full ov and reml alg 0"
gcta64 --mlma --bfile ${tissue}_GRM_input --grm ${tissue}_geno_grm --autosome-num 29 --reml-alg 1 \
--pheno phenotype.txt --mpheno ${pheno_id}  --qcovar cvrt.txt --reml-maxit 2000 --reml-bendV  \
--thread-num 10 --out output_mlm/${gene_list[i]}.assoc
if [ $? -gt 0 ]; then
echo "Failed attempt of gcta gwas of ${gene_list[i]} with reml alg 1"
gcta64 --mlma --bfile ${tissue}_GRM_input --grm ${tissue}_geno_grm --autosome-num 29 --reml-alg 2 \
--pheno phenotype.txt --mpheno ${pheno_id}  --qcovar cvrt.txt --reml-maxit 2000 --reml-bendV  \
--thread-num 10 --out output_mlm/${gene_list[i]}.assoc
if [ $? -gt 0 ]; then
echo "Failed attempt of gcta gwas of ${gene_list[i]} with reml alg 2"
fi
fi
fi

fi
awk -v gene=${gene_list[i]} '{if(NR==1)print "Gene""\t"$0; else if($9<0.00001)print gene"\t"$0}' output_mlm/${gene_list[i]}.assoc.mlma > output_mlm/${gene_list[i]}.assoc.mlma.tmp
rm output_mlm/${gene_list[i]}.assoc.mlma
done



