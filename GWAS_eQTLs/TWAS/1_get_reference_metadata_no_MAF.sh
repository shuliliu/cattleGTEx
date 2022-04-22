#!/bin/bash


##This script follows procedures in https://github.com/hakyimlab/summary-gwas-imputation 


REPO="~/bin/summary-gwas-imputation/src"
DATA="~/RNA-seq/run7_yak_indicus"

#module load bcftools

#chr=`for chr in $(seq 1 29); do echo ${chr}; done | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

####variant selection/Variant metadata
#rm $DATA/parquet_run7/Chr${chr}_maf0.01_monoallelic_variants.txt.gz
python $REPO/get_reference_metadata.py \
-genotype $DATA/parquet_run7/Chr${chr}-Run7-TAU-Beagle-toDistribute.txt.gz \
-annotation $DATA/parquet_run7/Chr${chr}-Run7-TAU-Beagle-toDistribute.annot.txt.gz \
-filter TOP_CHR_POS_BY_FREQ \
-rsid_column rsid \
-output $DATA/parquet_run7/Chr${chr}_monoallelic_variants.txt.gz

#Arrange the metadata file
zcat Chr1_monoallelic_variants.txt.gz | sed -n '1p' > variant_metadata_no_MAF.txt
for chr in $(seq 1 29)
do
gzip -cd Chr${chr}_monoallelic_variants.txt.gz | sed '1d' >> variant_metadata_no_MAF.txt
done
gzip -c variant_metadata_no_MAF.txt > variant_metadata_no_MAF.txt.gz

python $REPO/model_training_genotype_to_parquet.py \
-input_genotype_file $DATA/parquet_run7/Run7-TAU-Beagle-toDistribute.txt.gz \
-snp_annotation_file $DATA/parquet_run7/variant_metadata_no_MAF.txt.gz METADATA \
-parsimony 9 \
--impute_to_mean \
--only_in_key \
--split_by_chromosome \
-rsid_column rsid \
-output_prefix $DATA/parquet_run7/ARS_UCD1.2__no_MAF_monoallelic_variants

