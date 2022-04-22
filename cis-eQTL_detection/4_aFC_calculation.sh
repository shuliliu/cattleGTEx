#!/bin/bash

aFC="~/bin/aFC-master/aFC.py"

tissue="Adipose"
bcftools view -q 0.01:minor -v snps -c 4:minor ${tissue}.vcf > ./${tissue}/${tissue}.filtered.vcf
cd ${tissue}

Rscript ~/prepare_afc_phenotype.r ${tissue}

sort -k1,1 -k2,2n ${tissue}.phenotype.afc.bed | bgzip -cf > ${tissue}.phenotype.afc.bed.gz
rm ${tissue}.phenotype.afc.bed


sed '1d' ${tissue}.permutation.storey.txt | awk -F "_|\t" '{print $1"\t"$6"_"$7"_"$8"_"$9"\t"$6"\t"$7}' | sed '1 i pid\tsid\tsid_chr\tsid_pos' > ${tissue}.permutation.sig.txt

${aFC} --vcf ${tissue}.filtered.vcf --pheno ${tissue}.phenotype.afc.bed.gz --cov ${tissue}.covariance.txt.gz --gtl ${tissue}.permutation.sig.txt --output ${tissue}.permutation.sig.afc.txt

