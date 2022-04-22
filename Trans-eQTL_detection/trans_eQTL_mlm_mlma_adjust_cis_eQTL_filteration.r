tissue=commandArgs(T)[1]
ngene=as.numeric(commandArgs(T)[2])
nsnp=as.numeric(commandArgs(T)[3])

All_eQTL<-read.table(paste0("All_association.",tissue,".mlma_corrected.eQTL.txt"),header=T)

####Get the FDR.
All_eQTL<-All_eQTL[order(All_eQTL$P),]
All_eQTL$rank<-seq(1,nrow(All_eQTL))


All_eQTL$FDR<-All_eQTL$P*ngene*nsnp/(All_eQTL$rank)
All_eQTL_sig<-All_eQTL[All_eQTL$FDR<0.05,]

####Remove the cis-eQTL
used_genes<-read.table("~/ARS-UCD1.2/ARS1.2.genebody.protein_coding-lncRNA.bed",sep="\t")
All_eQTL_sig$gene_chr<-used_genes$V1[match(All_eQTL_sig$Gene,used_genes$V7)]
Trans_eQTL<-All_eQTL_sig[All_eQTL_sig$CHR!=All_eQTL_sig$gene_chr,]


####################################################################
###Cross Mappability remove
library("data.table")
cross_mappability<-read.table("~/cross_mappability3/cross_mappability.1Mb.genebody.txt")
colnames(cross_mappability)<-c("gene1","gene2","gene2_chr","gene2_1mbup","gene2_1mbdown")
cross_mappability<-cross_mappability[,c("gene1","gene2_chr","gene2_1mbup","gene2_1mbdown")]
step2<-Trans_eQTL[,c("SNP","Gene")]
filter<-merge(step2,cross_mappability,by.x="Gene",by.y="gene1")
filter<-as.data.frame(filter)
filter0<-as.data.frame(setDT(filter)[,c("gene1_chr","gene1_snp"):=tstrsplit(SNP,"_",keep=c(1,2))])
filter_genes<-filter[gene1_chr==gene2_chr & gene1_snp>=gene2_1mbup & gene1_snp<=gene2_1mbdown,]
gene_snp_pair<-filter_genes[,c("SNP","Gene")]

Trans_eQTL_filter<-Trans_eQTL[!paste(Trans_eQTL$SNP,Trans_eQTL$Gene)%in%paste(gene_snp_pair$SNP, gene_snp_pair$Gene),]
#######################################################################
length(unique(Trans_eQTL_filter$Gene))
length(unique(Trans_eQTL_filter$SNP))

write.table(Trans_eQTL_filter, "Trans_eQTL_FDR_0.05_filter_corrected_cis_eQTLs.txt",col.names=T, row.names=F, quote=F, sep="\t")


