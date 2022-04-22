#!/bin/bash

tissue="Adipose"
cd ${tissue}
mkdir dap_output
#We directly take the GTEx distributed single-SNP testing output to estimate the fine-mapping priors. In particular, this computation takes into account of SNP distance to transcription start site (TSS). The command to run is
torus="/home/shuli.liu/bin/dap-1.0.0/torus-master/src/torus"
cd ${tissue}

cd ../
${tissue}.nominals.txt.gz | awk '{print $1,$2,$3,$4,$5,0.05}OFS="\t"'  | gzip -c > ${tissue}.nominals.2rd2.txt.gz 
${torus} -d ${tissue}.nominals.2rd2.txt.gz --fastqtl -dump_prior ~/prior
rm ${tissue}.nominals.2rd2.txt.gz
