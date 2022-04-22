setwd("~/ASE/Adipose")



##### Only for top significant sQTLs per gene

aFC_slope<-read.table("Adipose_slope.ASE_intron.correlation_top_sQTL.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor.test(aFC_slope$afc,aFC_slope$slope,method="spearman") # 
corr
pdf(paste0("Adipose.aFC_slope_top_eQTL_intron.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "sQTL: Adipose, rho = -0.22, p = 4.09e-12")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()


##### Only for top significant sQTLs per gene

aFC_slope<-read.table("Muscle_slope.ASE_intron.correlation_top_sQTL.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor.test(aFC_slope$afc,aFC_slope$slope,method="spearman") # 
corr
pdf(paste0("Muscle.aFC_slope_top_eQTL_intron.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "sQTL: Muscle, rho = -0.229, p = 5.98e-27")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()



##### Only for top significant sQTLs per gene

aFC_slope<-read.table("Liver_slope.ASE_intron.correlation_top_sQTL.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor.test(aFC_slope$afc,aFC_slope$slope,method="spearman") # -0.357
corr
pdf(paste0("Liver.aFC_slope_top_eQTL_intron.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "sQTL: Liver, rho = -0.357, p = 0")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()