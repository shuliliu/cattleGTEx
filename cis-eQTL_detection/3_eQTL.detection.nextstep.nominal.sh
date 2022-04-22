#!/bin/bash

tissue="Adipose"

for j in $(seq 1 60); 
do
fastQTL.static --vcf ${tissue}.filtered.vcf.gz --bed ${tissue}.phenotype.bed.gz --cov ${tissue}.covariance.txt.gz  --normal --out ./Nominal/${tissue}.nominals.chunk${j}.txt.gz --chunk $j 60& 
done
wait

cd ./Nominal
zcat ${tissue}.nominals.chunk*.txt.gz | gzip -c > ${tissue}.nominals.txt.gz
rm ${tissue}.nominals.chunk*.txt.gz
Rscript ~/eQTL.p-value.nominal.correction_basePermutation.r ../Permutation/${tissue}.permutations.txt.gz 0.05 ${tissue}.nominals.txt.gz ${tissue}.nominals.sig.txt
