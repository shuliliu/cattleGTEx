#########################################################
#############Tissue-specific splicing#############
#########################################################
##clean up R environment
rm(list = ls())

##load required libraries
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("foreach")) {
  install.packages("foreach", dependencies = TRUE)
  library(foreach)
}

if (!require("doParallel")) {
  install.packages("doParallel", dependencies = TRUE)
  library(doParallel)
}




####################################################
###Gene expression(TPM) among 111 tissues/cell types in cattle

##import gene expression data matrix and sample information matrix

PSI0<-read.table("0-All-tissue_new_perind.counts.gz.phen_chrAll_sd",header=T) #column, samples; row, introns
PSI0<-PSI[,-(ncol(PSI0))]

rownames(PSI0)<-paste(PSI0[,c(1:4)],sep="-")
PSI0<-PSI0[,-c(1:4)]
head(PSI0[,1:5]);dim(PSI0)

info<-read.table("data_info_tissue_breed.txt",header=T,sep="\t")
head(info);dim(info); str(info)
info$Tissue_sample<-paste(info$Tissue_class,info$Sample,sep="_")

PSI0=t(PSI0) #let samples to be row; let introns to be column
##Only keep the matched rows in the intron PSI matrix
PSI<-subset(PSI0,rownames(PSI0)%in%info$Tissue_Sample)
dim(PSI)

sub_info<-subset(info,info$Tissue_Sample%in%rownames(PSI))
dim(sub_info)
sub_info<-sub_info[order(sub_info$Sample),]
PSI<-PSI[order(rownames(PSI)),]

save(PSI,sub_info,file="PSI-sub_info.RData")

##Summary of tissue and cell types
table(sub_info$Tissue_class)
length(unique(sub_info$Tissue_class))

##############################################################
############To identify the tissue specific genes#############

nsample=nrow(PSI)
ngene=ncol(PSI)

##step 1. write down the function for calculate the t value.
#T-statistic
#y=mean+x+e; X=1 (tested tissue); =-1 (tissues out of the tested system)
#FUNCTION
t_statistic<-function(y,x,Study,Breed,Sex,Age){
  linearMod<-lm(y~x+factor(Breed)+factor(Study)+factor(Sex)+Age)
  modelSummary <- summary(linearMod)  # capture model summary as an object
  modelCoeffs <- modelSummary$coefficients  # model coefficients
  xbeta.estimate <- modelCoeffs[2, "Estimate"]  # get beta estimate for speed
  std.error <- modelCoeffs[2, "Std. Error"]  # get std.error for speed
  t_value <- beta.estimate/std.error
  return(t_value)
}

p_statistic<-function(y,x,Study,Breed,Sex,Age){
  linearMod<-lm(y~x+factor(Breed)+factor(Study)+factor(Sex)+Age)
  modelSummary <- summary(linearMod)  # capture model summary as an object
  modelCoeffs <- modelSummary$coefficients  # model coefficients
  p_value<-modelCoeffs[2,"Pr(>|t|)"] # calc t statistic
  return(p_value)
}


##Step 2. Construct a matrix using loop to represent the existence of samples
ntissue=unique(sub_info$Tissue_class)
Matrix<-matrix(-1,nsample, length(ntissue))

j=1
for (i in ntissue){
  k=as.character(unique(sub_info[sub_info$Tissue_class==i,]$Tissue_categories))
  Matrix[sub_info$Tissue_class==i,j]<-1
  Matrix[which(sub_info$Tissue_categories==k & sub_info$Tissue_class!=i),j]<-NA
  j=j+1
}

colnames(Matrix)<-ntissue
rownames(Matrix)<-rownames(PSI)

dim(Matrix); class(Matrix); head(Matrix[,0:5])

##Step3. construct the the list of result t-value-matrix and p-value-matrix

t_value_matrix<-matrix(NA,ngene,length(ntissue))
colnames(t_value_matrix)<-ntissue
rownames(t_value_matrix)<-colnames(PSI)

dim(t_value_matrix); class(t_value_matrix); head(t_value_matrix[,0:5])

p_value_matrix<-matrix(NA,ngene,length(ntissue))
colnames(p_value_matrix)<-ntissue
rownames(p_value_matrix)<-colnames(PSI)


dim(p_value_matrix); class(p_value_matrix); head(p_value_matrix[,0:5])

##calculate the t-statistics
Sex=as.character(sub_info$Sex)
Sex[Sex=="Unknown"]<-"NA"

str(Sex)

Breed=as.character(sub_info$breed_class)
Breed[Breed=="Unknown"]<-"NA"

Age=sub_info$Age

Study=sub_info$Project

no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)


####
system.time({
for(j in ntissue){
t_value_matrix[,j]<-foreach (i=colnames(PSI),.combine=c)%dopar%{
t_statistic(PSI[,i],Matrix[,j],Study,Breed,Sex,Age) #gene expression, matrix, Study,Breed,Sex,Age
}

}
})
###

system.time({
  for(j in ntissue){
    p_value_matrix[,j]<-foreach (i=colnames(PSI),.combine=c)%dopar%{
      p_statistic(PSI[,i],Matrix[,j],Study,Breed,Sex,Age) #gene expression, matrix, Study,Breed,Sex,Age
    }
    
  }
})


write.csv(t_value_matrix,"t_statistic_114_tissues.matrix.splicing.csv")
write.csv(p_value_matrix,"p_statistic_114_tissues.matrix.splicing.csv")
stopImplicitCluster()




