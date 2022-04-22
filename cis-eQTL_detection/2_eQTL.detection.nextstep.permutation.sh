#!/bin/bash

tissue="Adipose"

for j in $(seq 1 60); 
do
fastQTL.static --vcf ${tissue}.filtered.vcf.gz --bed ${tissue}.phenotype.bed.gz --cov ${tissue}.covariance.txt.gz  --permute 1000 10000 --normal --out ./Permutation/${tissue}.permutations.chunk${j}.txt.gz --chunk $j 60& 
done
wait

cd ./Permutation
zcat ${tissue}.permutations.chunk*.txt.gz | gzip -c > ${tissue}.permutations.txt.gz
rm ${tissue}.permutations.chunk*.txt.gz
Rscript ~/eQTL.p-value.correction.r ${tissue} --no-save

