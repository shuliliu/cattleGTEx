
#######################################
##Plot eQTL
#input

gene=commandArgs(T)[1]
gene_name=commandArgs(T)[2]
tissue=commandArgs(T)[3]
trait=commandArgs(T)[4]
target_pos=as.numeric(commandArgs(T)[5])
target_snp=commandArgs(T)[6]
rcp=as.numeric(commandArgs(T)[7])
library(ggplot2)
library(ggrepel)
library(cowplot)

#format gwas
gwas<-read.table(paste0(trait, ".", tissue,".", gene,".gwas_r2.txt"))
colnames(gwas)<-c("chr","BP","logp","r2")
chr=unique(gwas$chr)

gwas$col[gwas$r2<=0.2]<-0.2
gwas$col[gwas$r2>0.2 & gwas$r2<=0.4]<-0.4
gwas$col[gwas$r2>0.4 & gwas$r2<=0.6]<-0.6
gwas$col[gwas$r2>0.6 & gwas$r2<=0.8]<-0.8
gwas$col[gwas$r2>0.8 & gwas$r2<=1]<-1

cols<-c("0.2"="blue4","0.4"="skyblue", "0.6"= "darkgreen","0.8"="orange","1"="red")  

gwas$target_snp<-target_snp

#format eqtl
eqtl<-read.table(paste0(trait, ".", tissue,".",gene,".eqtl_r2.txt"))
colnames(eqtl)<-c("chr","BP","logp","r2")

eqtl$col[eqtl$r2<=0.2]<-0.2
eqtl$col[eqtl$r2>0.2 & eqtl$r2<=0.4]<-0.4
eqtl$col[eqtl$r2>0.4 & eqtl$r2<=0.6]<-0.6
eqtl$col[eqtl$r2>0.6 & eqtl$r2<=0.8]<-0.8
eqtl$col[eqtl$r2>0.8 & eqtl$r2<=1]<-1

cols<-c("0.2"="blue4","0.4"="skyblue", "0.6"= "darkgreen","0.8"="orange","1"="red") 
eqtl$target_snp<-target_snp

##get the xlim 

start=min(min(eqtl$BP),min(gwas$BP))
end=max(max(eqtl$BP),max(gwas$BP))

##get the overlap between eqtl and gwas.

overlap<-merge(gwas,eqtl,by="BP")

rcp<-round(rcp,2)




#############################################################################
#plot combined plots.
subset=gwas[gwas$BP==target_pos,]
gwas_plot<-ggplot(gwas,aes(BP/1000000,logp,fill=factor(col)),color="black")+geom_point(shape=21,size=2)+
  geom_point(data=subset,aes(BP/1000000,logp),shape=23,size=4,fill="purple")+
  xlab("")+ylab(expression(GWAS-log(italic(P))))+
  theme_classic()+
  ggtitle(trait)+
  theme(        legend.position = "none",
                plot.title=element_text(size=16,color="black"),
                axis.title.y=element_text(size=16,color="black"),
                axis.title.x=element_text(size=16,color="black"),
                axis.text.y=element_text(size=16,color="black",hjust=0.5,vjust=0.5),
                axis.text.x=element_blank())+
  scale_fill_manual(values=cols)+
  #  geom_text_repel(
   # data = eqtl[eqtl$BP==target_pos,],
   # aes(label = target_snp),
  #  size = 3,
  #  box.padding = unit(0.35, "lines"),
  #  point.padding = unit(0.3, "lines")
  #)+
  #scale_y_continuous(breaks=c(0,1,2,3,4,5))+
  theme(plot.margin=unit(c(2,1.4,0.4,0.4),"cm"))#top right botton left

subset=eqtl[eqtl$BP==target_pos,]
xlabel=paste(chr,"position (Mb)")
eqtl_plot<-ggplot(eqtl,aes(BP/1000000,logp,fill=factor(col)),color="black")+geom_point(shape=21,size=2)+
  geom_point(data=subset,aes(BP/1000000,logp),shape=23,size=4,fill="purple")+
  xlab(xlabel)+ylab(expression(eQTL-log(italic(P))))+
  theme_classic()+
  ggtitle(paste0(tissue,": ",gene_name))+
  theme(        legend.position = "none",
                plot.title=element_text(size=16,color="black"),
                axis.title.y=element_text(size=16,color="black"),
                axis.title.x=element_text(size=16,color="black"),
                axis.text.y=element_text(size=16,color="black",hjust=0.5,vjust=0.5),
                axis.text.x=element_text(size=16,color="black",hjust=0.5,vjust=0.1))+
  scale_fill_manual(values=cols)+
  #geom_text_repel(
  #  data = eqtl[eqtl$BP==target_pos,],
  #  aes(label = target_snp),
  #  size = 3,
   # box.padding = unit(0.35, "lines"),
  #  point.padding = unit(0.3, "lines")
  #)+
  #scale_y_continuous(breaks=c(0,1,2,3,4,5))+
  theme(plot.margin=unit(c(0.2,1.4,0.4,0.4),"cm"))#top right botton left

subset=overlap[overlap$BP==target_pos,]
res=cor.test(overlap$logp.y, overlap$logp.x)
p=round(res$p.value,3)
r=round(res$estimate,2)

overlap_plot<-ggplot(overlap,aes(logp.y,logp.x),color="black")+geom_point(aes(fill=factor(col.x)),shape=21,size=2)+ #x eqtl; y gwas
  geom_point(data=subset,aes(logp.y,logp.x),shape=23,size=4,fill="purple")+
  xlab(expression(eQTL-log(italic(P))))+ylab(expression(GWAS-log(italic(P))))+

  theme_classic()+
  ggtitle(paste0("rcp = ", rcp))+
  theme(        plot.title=element_text(size=16,color="black"),
                axis.title.y=element_text(size=16,color="black"),
                axis.title.x=element_text(size=16,color="black"),
                axis.text.y=element_text(size=16,color="black",hjust=0.5,vjust=0.5),
                axis.text.x=element_text(size=16,color="black",hjust=0.5,vjust=0.1),
                legend.position="none")+
  scale_fill_manual(values=cols)+
  guides(fill=guide_legend(override.aes=list(shape=15,col=cols,size=4), title=bquote(~r^2)))+
   # geom_smooth(method=lm,se=FALSE)+
  #geom_text_repel(
  #  data = overlap[overlap$BP==target_pos,],
  #  aes(label = target_snp.x),
  #  size = 3,
   # box.padding = unit(0.35, "lines"),
   # point.padding = unit(0.3, "lines")
  #)+
  theme(plot.margin=unit(c(2,3,0.4,0.4),"cm"))
  #top right botton left
legend_box = data.frame(x = 0.8, y = seq(0.4, 0.28, -0.03))

overlap_plot=ggdraw(overlap_plot)+
             geom_rect(data = legend_box,
                      aes(xmin = x, xmax = x + 0.03, ymin = y, ymax = y + 0.03),
                      color = "black",
                      fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))) +
            draw_label("0.8", x = legend_box$x[1] + 0.03, y = legend_box$y[1], hjust = -0.3, size = 14) +
            draw_label("0.6", x = legend_box$x[2] + 0.03, y = legend_box$y[2], hjust = -0.3, size = 14) +
            draw_label("0.4", x = legend_box$x[3] + 0.03, y = legend_box$y[3], hjust = -0.3, size = 14) +
            draw_label("0.2", x = legend_box$x[4] + 0.03, y = legend_box$y[4], hjust = -0.3, size = 14) +
            draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.03, y = legend_box$y[1], vjust = -2, size = 14)

#####################################
#combine multiple plot into one plot.

library(ggpubr)

pdf(paste0("./Plots_with_cor/",trait, "_",tissue, "_",gene,"_",gene_name,".",target_snp,".plot.pdf"),width=10,height=6)

ggarrange(ggarrange(gwas_plot, eqtl_plot, nrow = 2, labels = c("A", "B"),font.label=list(color="black",face = "bold", size=18)),
           ggarrange(overlap_plot, labels = "C",font.label=list(color="black",face = "bold", size=18)), # Second row with box and dot plots
          ncol = 2                                      # Labels of the scatter plot
          ) 

dev.off()