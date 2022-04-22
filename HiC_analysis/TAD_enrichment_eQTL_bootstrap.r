library(data.table)
tissue=commandArgs(T)[1]
eQTL_dir="~/eQTL/vcf_tissues/"
out_dir="~/Hi-C/TAD_enrich_eQTL/"
gene_permut<-read.table(paste0(eQTL_dir,tissue,"/",tissue,".permutations.2rd.txt"))
gene_list<-gene_permut$V1 

gene_nominals<-read.table(paste0(eQTL_dir,tissue,"/",tissue,".nominals.2rd.sig.txt"))
gene_region<-read.table("~/ARS-UCD1.2/ARS1.2.genebody.bed")
colnames(gene_nominals)<-c("gene","nvar","dis","pval","slope","threshold")
colnames(gene_region)<-c("chr","start","end","gene","score","strand","cate")

###get the gene position and TSS

nominals_region<-merge(gene_nominals,gene_region,by="gene")


####get the subset of the gene-SNP pairs
if(nrow(nominals_region)>10000){
nominals_region_sub<-nominals_region[sample(nrow(nominals_region),10000),]}else{
nominals_region_sub<-nominals_region}

nominals_region_sub<-as.data.frame(setDT(nominals_region_sub)[,c("snp_loci","snp_ref","snp_alt"):=tstrsplit(nvar,"_",keep=c(2,3,4))])

nominals_region_sub$snp_loci<-as.numeric(nominals_region_sub$snp_loci)

nominals_region_sub$distance[nominals_region_sub$strand=="-"]<-nominals_region_sub$snp_loci[nominals_region_sub$strand=="-"]-nominals_region_sub$end[nominals_region_sub$strand=="-"]
nominals_region_sub$distance[nominals_region_sub$strand=="+"]<-nominals_region_sub$start[nominals_region_sub$strand=="+"]-nominals_region_sub$snp_loci[nominals_region_sub$strand=="+"] ##upstream of the TSSï¼? positive; downstream of the TSS: negative


###get the random gene-SNP pairs compatible to nominals_region_sub



random_pair<-function(num,gene_list,nominals_region_sub,gene_region){
for (i in 1:nrow(nominals_region_sub)){
    gene_random<-NULL
    nominals_region_sub$gene_random[i]<-as.character(gene_list[sample(length(gene_list),1)])

}

###get the random gene regions

nominals_region_random<-merge(nominals_region_sub,gene_region,by.x="gene_random",by.y="gene")
nominals_region_random$snp_loci.y[nominals_region_random$strand.y=="-"]<-nominals_region_random$end.y[nominals_region_random$strand.y=="-"]+nominals_region_random$distance[nominals_region_random$strand.y=="-"]
nominals_region_random$snp_loci.y[nominals_region_random$strand.y=="+"]<-nominals_region_random$start.y[nominals_region_random$strand.y=="+"]-nominals_region_random$distance[nominals_region_random$strand.y=="+"]

###
write_region_random<-nominals_region_random[,c("chr.y","start.y","end.y","gene_random","snp_loci.y")]

##write the bed format tables

write.table(write_region_random, paste0(out_dir,tissue,"/",tissue,".eGene.random.txt.tmp"),col.names=F,row.names=F,quote=F,sep="\t")

##run the overlap between TAD and eGene regions, to see which eGenes overlap


system(paste0("bedtools intersect -a ",out_dir,tissue,"/",tissue,".eGene.random.txt.tmp -b ~/Hi-C/hicMatrix/Lung_10kb_corrected_domains.bed -wa -wb > ",out_dir,tissue,"/",tissue,".eGene.random.TAD.txt.tmp"))

###to see whether the snps within the TAD regions
setwd(paste0(out_dir,tissue))
#system(paste0("head ",tissue,".eGene.random.TAD.txt.tmp"))


egene_random<-read.table(paste0(tissue,".eGene.random.TAD.txt.tmp"))
colnames(egene_random)<-c("chr","start","end","gene","snp_loci","TAD_chr","TAD_start","TAD_end","TAD_id")


egene_random$within<-"no"
egene_random$within[egene_random$snp_loci>=egene_random$TAD_start&egene_random$snp_loci<=egene_random$TAD_end]<-"yes"

summary(as.factor(egene_random$within))
egene_random_sum<-as.data.frame(summary(as.factor(egene_random$within)))
colnames(egene_random_sum)<-num
return(t(egene_random_sum))
}

eGene_random<-NULL
for (i in seq(1,5000)){
    tmp<-random_pair(i,gene_list,nominals_region_sub,gene_region)
eGene_random<-rbind(eGene_random,tmp)
}

write.table(eGene_random,paste0(out_dir,tissue,"/",tissue,".eGene_random_bootstrap.txt"),col.names=T,row.names=T,quote=F,sep="\t")

