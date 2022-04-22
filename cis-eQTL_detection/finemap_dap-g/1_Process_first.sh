#!/bin/bash

tissue="Adipose"
bcftools view -q 0.01:minor -v snps -c 4:minor ${tissue}.vcf > ./${tissue}/${tissue}.filtered.vcf
cd ${tissue}
bgzip -cf ${tissue}.filtered.vcf > ${tissue}.filtered.vcf.gz
tabix -p vcf ${tissue}.filtered.vcf.gz

##Converting to sbams format
#Run process script
process="~/dap-g_eQTL/process.pl"

perl ${process} -e ${tissue}.phenotype.bed.gz -g ${tissue}.filtered.vcf.gz -c ${tissue}.covariance.txt -t ${tissue}
