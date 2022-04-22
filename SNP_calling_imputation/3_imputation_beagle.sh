#!/bin/bash

Beagle="~/beagle.21Sep19.ec3.jar"
run7="~/run7/genotypes"

chrom=chr1

java -Xmx640g -Djava.io.tmpdir=$TMPDIR -XX:ParallelGCThreads=40 -jar ${Beagle} \
gt=cGTEx.${chrom}.vcf.gz \
ref=${run7}/Chr${chrom}-Run7-TAU-Beagle-toDistribute.vcf.gz \
out=../After_imputation/cGTEx.${chrom}-beagle.vcf.gz \
impute=true seed=9823 nthreads=40