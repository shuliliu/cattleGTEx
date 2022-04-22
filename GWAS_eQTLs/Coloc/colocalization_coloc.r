#!/bin/bash/Rscript

###input is eqtl file; gwas_file; tissue; trait; output_file;
args<-commandArgs(trailingOnly = T)
eQTL_file<-args[1]
GWAS_file<-args[2]
tissue=args[3]
trait=args[4]
output_dir=args[5]

library(coloc)
library(dplyr)
library(data.table)
if (!require("foreach")) {
  install.packages("foreach", dependencies = TRUE)
  library(foreach)
}

if (!require("doParallel")) {
  install.packages("doParallel", dependencies = TRUE)
  library(doParallel)
}

no_cores=4
cl<-makeCluster(no_cores)
registerDoParallel(cl)


gwas<-fread(GWAS_file,header=T)
colnames(gwas)<-c("SNP","frequency","pvalue","sample_size")
gwas_N=unique(gwas$sample_size)



  eqtl<-read.table(eQTL_file,header=T,stringsAsFactors=F)
  colnames(eqtl)<-c("Gene","SNP","Distance","P_value","Beta","Cutoff","N_sample","Freq_ref","Freq_alt")
  RES <- c("Tissue","Trait","Gene","nsnps","PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")
  eqtl$MAF=min(eqtl$Freq_ref,eqtl$Freq_alt)
  gene<-unique(eqtl$Gene)
  eqtl_N=unique(eqtl$N_sample)
  output_file=paste0(output_dir,"/",trait,".",tissue,".coloc.txt")


print(paste0("For ",tissue, " and ", trait, ", ", "Total genes: ", length(gene), " ; Begin process coloc......"))

system.time({
  #res_all<-foreach(j=gene,.combine=rbind)%dopar%{
    for (j in gene){
    library(coloc)
    print(j)
    qtl<-subset(eqtl,Gene==j)
  input <- merge(qtl, gwas, by="SNP", all=FALSE, suffixes=c("_eqtl","_gwas"))
  if(dim(input)[1]>0){
    result <- coloc.abf(dataset1=list(pvalues=input$P_value, type="quant", N=eqtl_N, MAF=input$MAF), dataset2=list(pvalues=input$pvalue, type="quant", N=gwas_N,MAF=input$frequency))
    res1 <- c(qtl$Gene[1])
    res2 <- as.data.frame(result$summary)
    res<-c(tissue,trait,res1,res2$`result$summary`)
    RES<-rbind(RES,res)
  }
}
})




write.table(RES,output_file,quote=F,col.names=F,row.names=F)
