#!/bin/bash



fastenloc="~/bin/fastenloc/src"
eQTL_dir="~/eQTL/vcf_tissues"
out_dir="~/colocalization/Cattle_GTEx_eQTL"

cd ${eQTL_dir}
tissue="Adipose"

perl ${fastenloc}/summarize_dap2enloc.pl -dir ${eQTL_dir}/${tissue}/dap_output -vcf ${eQTL_dir}/${tissue}/${tissue}.filtered.vcf.gz -tissue ${tissue} | gzip - > ${out_dir}/${tissue}.fastenloc.eqtl.annotation.vcf.gz
