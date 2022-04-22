#!/bin/bash

#### Work files
DIR_work=/home/lixj/data/analysis
DIR_download=/home/lixj/data
FILE_info_list=/home/lixj/data/pig_RNASeq.SRS_SRR.list

### Number of threads
NSLOTS=8
### Required in Step1: RNA_Seq_Pipeline_QC_Salmon_Align_Quant.sh
Reference_index_Salmon=/vloume01/pig_refgenome/pig_transcripts_index/
Reference_index_STAR=/vloume01/pig_refgenome/STAR_index/
Reference_GTF=/vloume01/pig_refgenome/Sus_scrofa.Sscrofa11.1.100.gtf

### Required in Step2: RNA_Seq_Pipeline_SNP_ASE_Splicing.sh
Reference_Fasta=/media/disk3/sijf/Buffalo_RNASeq/Buffalo_genome/GCF_003121395.1_ASM312139v1_genomic.fna
###???dbSNP=/media/disk3/sijf/Buffalo_RNASeq/Buffalo_genome/Buffalo_451_rmfilter_biallelic_concat_NC.snp.vcf.gz
interval=/vloume01/pig_refgenome/pig.autosome.list

### Softwares
TRIMMOMATICDIR=/home/lixj/software/Trimmomatic-0.39
PICARDDIR=/home/lixj/software
LeafCutterDIR=/home/lixj/software/leafcutter-master

# export PATH=/home/lixj/software/bwa-0.7.17:$PATH
# export PATH=/home/lixj/software/samtools-1.10:$PATH
# export PATH=/home/lixj/software/salmon-latest_linux_x86_64/bin:$PATH
# export PATH=/home/lixj/software/STAR-2.7.3a/bin/Linux_x86_64:$PATH
# export PATH=/home/lixj/software/stringtie-2.1.2.Linux_x86_64:$PATH
# export PATH=/home/lixj/software/subread-2.0.0-Linux-x86_64/bin:$PATH
# export PATH=/home/lixj/software/gatk-4.1.7.0:$PATH