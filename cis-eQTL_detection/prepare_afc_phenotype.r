moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

library(SNPRelate)
tissue=commandArgs(T)[1]
vcf.fn<-paste0(tissue,".filtered.vcf")
snpgdsVCF2GDS(vcf.fn,"ccm.gds",method="biallelic.only")
genofile<-openfn.gds("ccm.gds")
ccm_pca<-snpgdsPCA(genofile)

if(length(ccm_pca$sample.id) < 150 ) {
pca_genotype<-ccm_pca$eigenvect[,1:3]
colnames(pca_genotype)<-c("pc1","pc2","pc3")
}else if (length(ccm_pca$sample.id)  < 250){
pca_genotype<-ccm_pca$eigenvect[,1:5]
colnames(pca_genotype)<-c("pc1","pc2","pc3","pc4","pc5")
}else{
pca_genotype<-ccm_pca$eigenvect[,1:10]
colnames(pca_genotype)<-c("pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")
}

rownames(pca_genotype)<-ccm_pca$sample.id

#####expression######
#Phenotype matrix: 
#Chromosome ID [string]
#Start genomic position of the phenotype (e.g. TSS) [integer]
#End genomic position of the phenotype (e.g. TTS) [integer]
#Phenotype ID [string] Then additional columns give phenotype quantifications for all samples. 
library(peer)
library(preprocessCore)
library(RNOmni)
load("~/RNA-seq/all.LJA.8742samples.RData") 
expr=subset(TPM_tmp0,rownames(TPM_tmp0)%in%ccm_pca$sample.id)
expr_matrix00=as.data.frame(t(expr))
expr_matrix00<-expr_matrix00[,order(factor(colnames(expr_matrix00),levels=ccm_pca$sample.id))]
count_0.1<-rowSums(expr_matrix00>0.1)
#calculate the sample number:
nsamples=length(ccm_pca$sample.id)
#keep the genes with >0.1 tpm in >=20% samples.
expr_matrix00<-expr_matrix00[count_0.1>=(0.2*nsamples),] ##row is gene; column is sample
expr_matrix00$id<-row.names(expr_matrix00)
region_annot<-read.table("~/ARS-UCD1.2/ARS1.2.TSS.bed")
colnames(region_annot)<-c("#Chr","start","end","name","score","strand","ID")
expr_matrix0<-merge(region_annot,expr_matrix00,by.x="ID",by.y="id")
expr_matrix<-expr_matrix0[,-which(names(expr_matrix0) %in% c("name","score","strand"))]
expr_matrix<-expr_matrix[moveme(names(expr_matrix),"ID after end")]
##expr_matrix 
write.table(expr_matrix,paste0(tissue,".phenotype.afc.bed"),sep="\t",row.names=F,quote =FALSE)
