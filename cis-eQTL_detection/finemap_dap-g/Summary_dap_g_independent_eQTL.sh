#!/bin/bash

cd ~/eQTL/vcf_tissues
tissue="Adipose"
bcftools view -q 0.01:minor -v snps -c 4:minor ${tissue}.vcf > ./${tissue}/${tissue}.filtered.vcf
cd ${tissue}

bgzip -cf ${tissue}.filtered.vcf > ${tissue}.filtered.vcf.gz
tabix -p vcf ${tissue}.filtered.vcf.gz

##Converting to sbams format
#Run process script
process="/home/shuli.liu/commands/commands_for_RNA_seq_2/eQTL/dap-g_eQTL/process.pl"

perl ${process} -e ${tissue}.phenotype.bed.gz -g ${tissue}.filtered.vcf.gz -c ${tissue}.covariance.txt -t ${tissue}

#sh Adipose.assemble.cmd
##Estimate priors for fine-mapping
#We directly take the GTEx distributed single-SNP testing output to estimate the fine-mapping priors. In particular, this computation takes into account of SNP distance to transcription start site (TSS). The command to run is
torus="/home/shuli.liu/bin/dap-1.0.0/torus_src/torus"
dap="/home/shuli.liu/bin/dap-1.0.0/dap-master/dap_src/dap-g"
cd ${tissue}
mkdir prior
cd ../
${torus} -d ${tissue}.nominals.2rd.txt.gz -dump_prior ${tissue}/prior

##The command we use for GTEx analysis is
${dap} -d ./Adipose/ENSBTAG00000019404.sbams.dat -p ./Adipose/prior/ENSBTAG00000019404.prior -ld_control 0.5 --all -t 4 > ENSBTAG00000019404.dap.txt
#-d and -p point to the sbams and prior files. Option --all forces outputting information of all SNPs (not just noteworthy ones). Note that --all has now become a default option in DAP. -t 4 indicates the run will use 4 parallel threads if available. -ld_control 0.5 specifies the lowest LD threshold (R^2) to admit a SNP into a signal cluster.

