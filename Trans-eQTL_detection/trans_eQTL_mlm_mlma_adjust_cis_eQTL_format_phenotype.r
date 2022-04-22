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
finemap_dir="~/vcf_tissues/fine_map_eQTL/"

####Format the genotypes ###We don't need to format the genotype anymore

####Format the covariances
####I need to consider two situiations: 1. gene without finemapping result; 2. gene with finemapping result. it seems for each gene we should generate a fixed effect file. 
snp_file<-read.table(paste0(trans_orig,"/",tissue, "/",tissue, ".snp.txt"),header=T)
gene_list<-read.table("gene_list.txt")
finemap=read.table(paste0(finemap_dir,tissue,".finemap.eQTL.txt"))

cvrt<-read.table(paste0(trans_orig,"/",tissue, "/",tissue, ".cvrt.txt"),header=T)
rownames(cvrt)<-cvrt$id
cvrt<-as.data.frame(t(cvrt[,-1]))

####Make a directory storing the fixed effect files.

system("mkdir fixed")


if ("Species" %in% colnames(cvrt)){  ##check whether have fixed effect: species.

fixed_species<-as.data.frame(cvrt[1]) 
fixed_species$id<-rownames(cvrt)
fixed_species$id2<-rownames(cvrt)
fixed_species<-fixed_species[moveme(names(fixed_species),"id first")]
fixed_species<-fixed_species[moveme(names(fixed_species),"id2 first")]

for (gene in gene_list$V1){   ## for each gene, generate a fixed effect file.

if (gene%in%finemap$V2){
snps<-finemap[finemap$V2==gene,]$V1
genotype<-snp_file[snp_file$id%in%snps,]
rownames(genotype)<-genotype$id
genotype<-t(genotype[,-1])
fixed<-cbind(fixed_species, genotype)
}else{
fixed<-fixed_species
}
write.table(as.data.frame(fixed),paste0("./fixed/", gene,".fixed.txt"),col.names=F,row.names=F,quote=F,sep="\t")
}

}else{
for (gene in gene_list$V1){
if (gene%in%finemap$V2){
snps<-finemap[finemap$V2==gene,]$V1
genotype<-snp_file[snp_file$id%in%snps,]
rownames(genotype)<-genotype$id
genotype<-t(genotype[,-1])
genotype$id<-rownames(genotype)
    genotype$id2<-rownames(genotype)
genotype<-genotype[moveme(names(genotype),"id first")]
genotype<-genotype[moveme(names(genotype),"id2 first")]
write.table(as.data.frame(fixed),paste0("./fixed", gene,".fixed.txt"),col.names=F,row.names=F,quote=F,sep="\t")
}
}
}

