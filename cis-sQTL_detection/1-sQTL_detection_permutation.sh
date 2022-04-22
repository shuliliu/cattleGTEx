# Take Liver as an example

# Step 1. Converting bams to juncs
for i in $(seq 1 10)
do
	sh leafcutter/scripts/bam2junc.sh liver${i}.bam liver${i}.junc
done


# Step 2. Intron clustering
ls *.junc > liver_juncfiles.txt
python leafcutter/clustering/leafcutter_cluster.py -j liver_juncfiles.txt -m 50 -o liver -l 500000


# Step 3. Calculate the PCA 
python leafcutter/scripts/prepare_phenotype_table.py liver_perind.counts.gz -p 10 # note: -p 10 specific you want to calculate for sQTL to use as covariates
sh liver_perind.counts.gz_prepare.sh
head -6 liver_perind.counts.gz.PCs > liver_perind.counts.gz.PCs.PC5


# Step 4. sQTL detection
for j in $(seq 1 29)
do
	fastQTL.static --vcf liver.filtered.vcf.gz --bed liver_perind.counts.gz.qqnorm_chr${j}.gz --cov liver_perind.counts.gz.PCs.PC5 --permute 1000 10000 --normal --out liver_perind.permutation.chr${j} --chunk 1 1
done

# Step 5. Calculate the FDR using Benjamini & Hochberg method and storey-Tibshriani method.
Rscript FDR_correction.r




