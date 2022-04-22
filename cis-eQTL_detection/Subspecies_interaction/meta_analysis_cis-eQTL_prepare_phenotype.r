###function#########
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

tissue=commandArgs(T)[1]
gene=commandArgs(T)[2]
#######################################
#Prepare the genotypes
######################################
setwd(paste0("~/eQTL/vcf_tissues/",tissue))
system(paste0("sed 's/#//g' ",tissue,".phenotype.bed > phenotype.bed.tmp"))
phenotype<-read.table("phenotype.bed.tmp",header=T)
system("rm phenotype.bed.tmp")
phenotype_gene<-as.vector(t(phenotype[phenotype$ID==gene,-c(1:4)]))


###format the corrected phenotype FID FID genotype
co_pheno<-as.data.frame(phenotype_gene)
rownames(co_pheno)<-colnames(phenotype[,-c(1:4)])
colnames(co_pheno)<-gene
co_pheno$IID<-rownames(co_pheno)
co_pheno$FID<-rownames(co_pheno)
co_pheno<-co_pheno[moveme(names(co_pheno),"IID first")]
co_pheno<-co_pheno[moveme(names(co_pheno),"FID first")]

setwd(paste0("~/cis-eQTL_Muscle_test_SNP_breed/",tissue))
write.table(co_pheno, paste0(gene,".phenotype.txt.tmp"),col.names=T,row.names=F,quote=F,sep="\t")

###format the breed covariate.
#setwd(paste0("~/eQTL/vcf_tissues/",tissue))
#breed<-read.table(paste0(tissue,".covariance.txt"),header=T)
#rownames(breed)<-breed$id
#breed<-t(breed[,-1])
#breed<-as.data.frame(breed)
#breed$Species2[breed$Species=="Bos_taurus"]=0
#breed$Species2[breed$Species=="Bos_indicus"]=1
#breed$Species2[breed$Species=="Bos_indicus_x_Bos_taurus"]=2
#breed$Species2[breed$Species=="Bos_grunniens"]=3
#breed$Species2[breed$Species=="Bos_grunniens_x_Bos_taurus"]=4
#breed$Species2[breed$Species=="Bos_frontalis"]=5
#breed$Species2[breed$Species=="Buffalo"]=6
#breed<-breed[,!(colnames(breed)%in%c("Species"))]
#breed<-breed[moveme(names(breed),"Species2 first")]
#breed$IID<-rownames(breed)
#breed$FID<-rownames(breed)
#breed<-breed[moveme(names(breed),"IID first")]
#breed<-breed[moveme(names(breed),"FID first")]
#setwd(paste0("~/cis-eQTL_Muscle_test_SNP_breed/",tissue))
#write.table(breed, "covariate.txt",col.names=T,row.names=F,quote=F,sep="\t")


















