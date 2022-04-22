#!/bin/bash

tissue="Adipose"

###minor allele frequencies >=0.01; with the minor allele observed in at least 4 samples.
bcftools view -q 0.01:minor -v snps -c 4:minor ./vcf_tissues/${tissue}.vcf > ${tissue}.filtered.vcf
bgzip -cf ${tissue}.filtered.vcf > ${tissue}.filtered.vcf.gz
tabix -p vcf ${tissue}.filtered.vcf.gz

##prepare the phenotype and covariance
Rscript ~/prepare_for_eQTL.detection.r ${tissue} --no-save

##format them.
sort -k1,1 -k2,2n ${tissue}.phenotype.bed  | bgzip -cf > ${tissue}.phenotype.bed.gz
tabix -p bed ${tissue}.phenotype.bed.gz

bgzip -cf ${tissue}.covariance.txt > ${tissue}.covariance.txt.gz

