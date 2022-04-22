#!/bin/bash



module load r
module load bedtools

eQTL_dir="~/eQTL/vcf_tissues"
gene_region="~/ARS-UCD1.2/ARS1.2.genebody.bed"
TSS="~/ARS-UCD1.2/ARS1.2.TSS.bed"
aim_dir="~/Hi-C/TAD_enrich_eQTL"
cd ${eQTL_dir}
tissues=("Adipose" "Blood" "Embryo" "Hypothalamus" "Ileum" "Intramuscular_fat" "Jejunum" "Leukocyte" "Liver" "Lung" "Lymph_node"
"Macrophage" "Mammary" "Milk_cell" "Monocytes" "Muscle" "Ovary"
"Oviduct" "Pituitary" "Rumen" "Salivary_gland" "Skin_fibroblast"
"Testis" "Uterus")

index=$(($SLURM_ARRAY_TASK_ID-1))

tissue=${tissues[index]}

cd ${aim_dir}
mkdir ${tissue}
cd ${eQTL_dir}  

Rscript TAD_enrichment_eQTL_bootstrap.r ${tissue}


  
