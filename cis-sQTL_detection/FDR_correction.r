library(qvalue)
d <- read.table("liver_perind.permutation.chrAll.txt", header = F, stringsAsFactors=F)
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")

pdf(paste0(Tissue,"liver_perind.permutation.pdf"),width=6,height=6)
plot(d$ppval, d$bpval, xlab="Direct method", ylab="Beta approximation", main="Check plot")
abline(0, 1, col="red")
dev.off()

d$bh = p.adjust(d$bpval, method="fdr")
write.table(d[which(d$bh <= 0.05), ], "liver_perind.permutation.fdr.txt", quote = F, row.names=F, col.names = T, sep = "\t")
d$st = qvalue(d$bpval,lambda=0.5)$qvalues
write.table(d[which(d$st <= 0.05), ], "liver_perind.permutation.storey.txt", quote = F, row.names=F, col.names = T, sep = "\t")