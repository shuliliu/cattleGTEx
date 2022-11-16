setwd("/home/shuli.liu/uvm_mckay/WGBS_others/RNA-seq/TPM_matrix")
###Plot the cluster simiar to tsne to learn the hierarchical clustering of gene expression using our samples
library("dendextend")
library(RColorBrewer)
###
#from load part
load("../all.LJA.8742samples.RData")
####
info<-read.table("../data_info_tissue_breed_clean_reads_2M.txt",header=T,sep="\t")
TPM<-subset(TPM_tmp0,rownames(TPM_tmp0)%in%info$Sample)
dim(TPM)
sub_info<-subset(info,info$Sample%in%rownames(TPM))
dim(sub_info)
sub_info<-sub_info[order(sub_info$Sample),]
TPM<-TPM[order(rownames(TPM)),]

TPM_count<-colSums(TPM>0)
TPM0=TPM[,TPM_count!=0]
TPM0_count=TPM_count[TPM_count!=0]
TPM_loge_scale<-t(apply(log(TPM0+0.25), MARGIN = 1, FUN = scale)) #TPM_loge_scale: row is individuals, col is variable (gene expression)

###########
##calculate the Standard deviation

TPM_loge_scale_sd<-apply(TPM_loge_scale,MARGIN=2, sd)

##order the TPM_loge_scale based on sd (large to small)
TPM_loge_scale<-TPM_loge_scale[,order(TPM_loge_scale_sd,decreasing=T)]#Assign colors

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors<-col_vector[c(5,6,13,9,3,2,11,7,14,21,22,1,50)]
names(colors)=unique(sort(sub_info$Tissue_categories))
colors["Other"]="#E0E0E0"
colors["Pituitary"]="black"
colors["Kidney"]="#0F27F0"
sub_info$colors<-colors[match(sub_info$Tissue_categories,names(colors))]

#delete the others sample:
TPM_loge_scale<-TPM_loge_scale[sub_info$Tissue_categories!="Other",]
sub_info<-sub_info[sub_info$Tissue_categories!="Other",]

#calculate the dend 

dend <- as.dendrogram(hclust(dist(TPM_loge_scale)))  #row is individuals, col is the variable (gene expression)
save(dend,file="dend.RData")
  pdf("colored_bar.plot.TPM_loge_scale.pdf",width=12,height=5)
  par(mar = c(12,2,1,4))
  plot(dend,leaflab="none",axes=F,edge.color ="light grey")
  colored_bars(sub_info$colors, dend, y_shift=-4, y_scale=15, rowLabels = "Tissues",borders(fill = NA))
  dev.off()

