interaction<-read.table("ADDxSpecies_info.txt")
colnames(interaction)<-c("id","CHR","SNP","BP","A1", "TEST","NMISS", "BETA", "STAT","P")
interaction$FDR<-p.adjust(interaction$P,method="BH")
write.table(interaction, "ADDxSpecies_info.FDR.txt",col.names=T, row.names=F, quote=F, sep="\t")

