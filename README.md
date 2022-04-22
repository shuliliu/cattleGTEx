Cattle GTEx (v1) 
=================
# 1. Introduction
Characterization of genetic regulatory variants acting on the transcriptome of livestock is essential for interpreting the molecular mechanisms underlying traits of economic value and for increasing the rate of genetic gain through artificial selection. Here we build a **cattle Genotype-Tissue Expression atlas (cattleGTEx)** as part of the pilot phase of **Farm animal GTEx (FarmGTEx)** project for the research community based on publicly available **7,180** RNA-Seq samples. We describe the landscape of transcriptome across over 100 tissues and report hundreds of thousands of genetic associations with gene expression and alternative splicing for **24** major tissues. We evaluate the tissue-sharing patterns of these genetic regulatory effects, and functionally annotate them using multi-omics data. Finally, we link gene expression in different tissues to **43** economically important traits using both transcriptome-wide association and colocalization analyses to decipher the molecular regulatory mechanisms underpinning such agronomic traits in cattle. 

# 2. The main contents
There are mainly nine parts of analysis in our project. 

## Gene_expression_quantification
Includes gene expression quantification from RNA-seq.

## DNA_methylation
Includes the analysis of WGBS data.

## Tissue_specific_expression_splicing
Includes the detection of tissue specific expressed genes and spliced introns using the method illustrated in [Finucane et al. (2018)](https://www.nature.com/articles/s41588-018-0081-4).

## ASE
Includes the estimates of effect size (aFC) of requlatory variants using aggregated phASER haplotypic expression data using [phASER](https://github.com/secastel/phaser) and the correlation plot between top cis-eQTL slopes (from [fastQTL](https://github.com/francois-a/fastqtl)).

## SNP_calling_imputation
Includes the SNP calling from RNA-seq and imputation using [Beagle 5](https://faculty.washington.edu/browning/beagle/beagle.html).

## *cis*-eQTL_detection
Includes the eQTL detection using [fastQTL](https://github.com/francois-a/fastqtl), effect size calculation using [aFC](https://github.com/secastel/aFC) and finemapping analysis using [DAP-G](https://github.com/xqwen/dap)

## *Trans*-eQTL_detection
Includes the *trans*-eQTL detection by a mixed linear model using mlma from [GCTA](https://yanglab.westlake.edu.cn/software/gcta/#Overview) with and without cis-eQTL adjustments.

##  GWAS_eQTLs
Includes the TWAS analysis using S-PrediXcan and MultiXcan from [MetaXcan](https://github.com/hakyimlab/MetaXcan) and Colocalization analysis using [Coloc](https://github.com/chr1swallace/coloc) and [fastENLOC](https://github.com/xqwen/fastenloc).

## HiC_analysis
Includes the HiC data processes using [HiC-Pro(v2.11.4)](https://github.com/nservant/HiC-Pro) and estimates of significant intra-chromosome contacts using [FitHiC(v2.0.7)](https://github.com/ay-lab/fithic).



