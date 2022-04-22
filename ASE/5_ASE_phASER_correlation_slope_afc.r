setwd("~/ASE")

#####For all sig nominal eQTLs

aFC_slope<-read.table("Adipose_slope.ASE.correlation.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor.test(aFC_slope$afc,aFC_slope$slope,method="spearman") #0.6868973
corr
pdf(paste0("Adipose.aFC_slope.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "eQTL: Adipose, rho = 0.69, p = 0")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()

##### Only for top significant eQTLs per gene

aFC_slope<-read.table("Adipose_slope.ASE.correlation_top_eQTL.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor(aFC_slope$afc,aFC_slope$slope,method="spearman") # 0.7542782
corr
'''
data:  aFC_slope$afc and aFC_slope$slope
S = 1.5132e+12, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho
0.6868973
pvalue:0
'''

pdf(paste0("Adipose.aFC_slope_top_eQTL.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "eQTL: Adipose, rho = 0.75, p = 0")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()


####For Muscle
aFC_slope<-read.table("Muscle_slope.ASE.correlation_top_eQTL.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor.test(aFC_slope$afc,aFC_slope$slope,method="spearman") # 0.676
corr

'''
 Spearman rank correlation rho

data:  aFC_slope$afc and aFC_slope$slope
S = 95647648, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho
0.6760562
p.value: 2.110673e-162
'''

pdf(paste0("Muscle.aFC_slope_top_eQTL.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "eQTL: Muscle, rho = 0.676, p = 0")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()

####For Liver
aFC_slope<-read.table("Liver_slope.ASE.correlation_top_eQTL.txt")
colnames(aFC_slope)<-c("gene","var","afc","slope")
corr<-cor.test(aFC_slope$afc,aFC_slope$slope,method="spearman") # 0.735
corr
pdf(paste0("Liver.aFC_slope_top_eQTL.pdf"),width=4,height=4.5)
plot(aFC_slope$slope,aFC_slope$afc,xlab="FastQTL slope",ylab="ASE aFC",main = "eQTL: Liver, rho = 0.735, p = 0")
abline(lm(aFC_slope$afc~aFC_slope$slope),lty=2,col="dark grey")
dev.off()