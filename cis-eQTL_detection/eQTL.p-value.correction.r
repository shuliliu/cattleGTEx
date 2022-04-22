########checking that the experiment went well; controlling for multiple phenotypes besing tested#################
tissue=commandArgs(T)[1]

d=read.table(gzfile(paste0(tissue,".permutations.txt.gz")),hea=F,stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")
pdf(paste0(tissue,".check.p-values3.pdf"),width=6,height=6)
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
dev.off()

####Storey & Tibshirani correction
library(qvalue)
#d$st = qvalue(d$bpval,lamda=0.85)$qvalues
d$bh=p.adjust(d$bpval,method="BH")
write.table(d[which(d$bh <= 0.05), ], paste0(tissue,".permutations.storey.txt"), quote=F, row.names=F, col.names=T)
