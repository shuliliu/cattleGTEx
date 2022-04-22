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
trans_vcf="~/eQTL/trans"


####Format the genotypes
gene<-read.table(paste0(trans_vcf,"/",tissue, "/",tissue, ".gene.txt"),header=T)
gene_list<-gene$ID
rownames(gene)<-gene$ID
gene<-as.data.frame(t(gene[,-1]))
gene$id<-rownames(gene)
gene$id2<-rownames(gene)
gene<-gene[moveme(names(gene),"id first")]
gene<-gene[moveme(names(gene),"id2 first")]

write.table(as.data.frame(gene_list),"gene_list.txt",col.names=F, row.names=F, quote=F, sep="\t")
write.table(gene, "phenotype.txt",col.names=F,row.names=F,quote=F,sep="\t")

####Format the covariances
cvrt<-read.table(paste0(trans_orig,"/",tissue, "/",tissue, ".cvrt.txt"),header=T)
rownames(cvrt)<-cvrt$id
cvrt<-as.data.frame(t(cvrt[,-1]))
if ("Species" %in% colnames(cvrt)){
    fixed<-as.data.frame(cvrt[1])
    fixed$id<-rownames(cvrt)
    fixed$id2<-rownames(cvrt)
    fixed<-fixed[moveme(names(fixed),"id first")]
    fixed<-fixed[moveme(names(fixed),"id2 first")]
    write.table(as.data.frame(fixed),"fixed.txt",col.names=F,row.names=F,quote=F,sep="\t")
    cvrt<-cvrt[,-1]
}

cvrt$id<-rownames(cvrt)
cvrt$id2<-rownames(cvrt)
cvrt<-cvrt[moveme(names(cvrt),"id first")]
cvrt<-cvrt[moveme(names(cvrt),"id2 first")]
write.table(cvrt, "cvrt.txt",col.names=F,row.names=F,quote=F,sep="\t")
