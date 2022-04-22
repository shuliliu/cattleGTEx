#!/bin/bash

################filter the snps DR2>=0.8; minor allele frequency: 0.05; only keep the snps

chr=chr1

bcftools view -q 0.05:minor -v snps -i "DR2>=0.8" ${chr}.all-beagle.vcf.gz -Oz > ./filter_MAF_0.05_DR2_0.8/${chr}-beagle.filter.maf.vcf.gz

