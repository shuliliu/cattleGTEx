#!/bin/bash

fastQTL="~/bin/FastQTL/bin"
vcf="~/vcf_tissues"
bed="~/3-tissue-all-split"

Tissue="Liver"

cd ${Tissue}/nominal

for j in $(seq 1 29)
do
	for i in $(seq 1 300)
	do
		${fastQTL}/fastQTL.static --vcf ${vcf}/${Tissue}.vcf.gz --bed ${bed}/${Tissue}/${Tissue}.qqnorm_chr${j}.gz --cov ../${Tissue}_pca.txt --normal --out ${Tissue}.nominal.chr${j}.chunk${i}.txt --chunk $i 300
		cat ${Tissue}.nominal.chr${j}.chunk${i}.txt >> ${Tissue}.nominal.chr${j}.txt
	done
done

# note: when some issue like "0m Number of chunks (500) is greater than the number of phenotypes (394)" happens, if we get the expecting results, then we can ignore it, but if the results are empty,
# then we should set the number of chunks lower than the number listed in the (), say 394 here. For example, we can set 300.
