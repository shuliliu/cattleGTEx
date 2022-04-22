#########################################################
#############Tissue-specific gene expression#############
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
TPM_tmp0<-get(load("../all.LJA.8742samples.RData"))
head(TPM_tmp0[,1:5]);dim(TPM_tmp0)

info<-read.table("../data_info_tissue_breed.txt",header=T,sep="\t")
head(info);dim(info); str(info)

##Only keep the matched rows in the gene expression matrix
TPM<-subset(TPM_tmp0,rownames(TPM_tmp0)%in%info$Sample)
dim(TPM)

sub_info<-subset(info,info$Sample%in%rownames(TPM))
dim(sub_info)
sub_info<-sub_info[order(sub_info$Sample),]
TPM<-TPM[order(rownames(TPM)),]



##Summary of tissue and cell types
table(sub_info$Tissue_class)
length(unique(sub_info$Tissue_class))

##############################################################
############To identify the tissue specific genes#############


##remove not-expressed genes; center and scale all the columns, i.e., normalize across samples.
TPM_count<-colSums(TPM>0)
TPM0=TPM[,TPM_count!=0]
TPM0_count=TPM_count[TPM_count!=0]
TPM_log10_scale<-t(apply(log(TPM0+0.25), MARGIN = 1, FUN = scale))
colnames(TPM_log10_scale)<-colnames(TPM0)

dim(TPM_log10_scale); class(TPM_log10_scale); head(TPM_log10_scale[,0:5])


nsample=nrow(TPM_log10_scale)
ngene=ncol(TPM_log10_scale)

##step 1. write down the function for calculate the t value.
#T-statistic
#y=mean+x+e; X=1 (tested tissue); =-1 (tissues out of the tested system)
#FUNCTION
t_statistic_t<-function(y,x,Study,Breed,Sex,Age){
  linearMod<-lm(y~x+factor(Breed)+factor(Study)+factor(Sex)+Age)
  modelSummary <- summary(linearMod)  # capture model summary as an object
  modelCoeffs <- modelSummary$coefficients  # model coefficients
  beta.estimate <- modelCoeffs[2, "Estimate"]  # get beta estimate for speed
  std.error <- modelCoeffs[2, "Std. Error"]  # get std.error for speed
  t_value <- beta.estimate/std.error
  #p_value<-modelCoeffs[2,"Pr(>|t|)"] # calc t statistic
  return(t_value)
}

t_statistic_p<-function(y,x,Study,Breed,Sex,Age){
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
rownames(Matrix)<-rownames(TPM_log10_scale)

dim(Matrix); class(Matrix); head(Matrix[,0:5])

##Step3. construct the the list of result t-value-matrix and p-value-matrix

t_value_matrix<-matrix(NA,ngene,length(ntissue))
colnames(t_value_matrix)<-ntissue
rownames(t_value_matrix)<-colnames(TPM_log10_scale)

dim(t_value_matrix); class(t_value_matrix); head(t_value_matrix[,0:5])

p_value_matrix<-matrix(NA,ngene,length(ntissue))
colnames(p_value_matrix)<-ntissue
rownames(p_value_matrix)<-colnames(TPM_log10_scale)


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


####calculate the values
system.time({
  for(j in ntissue){
  t_value_matrix[,j]<-foreach (i=colnames(TPM_log10_scale),.combine=c)%dopar%{
    t_statistic_t(TPM_log10_scale[,i],Matrix[,j],Study,Breed,Sex,Age) #gene expression, matrix, Study,Breed,Sex,Age
    }
  
  }
})


system.time({
for(j in ntissue){
p_value_matrix[,j]<-foreach (i=colnames(TPM_log10_scale),.combine=c)%dopar%{
t_statistic_p(TPM_log10_scale[,i],Matrix[,j],Study,Breed,Sex,Age) #gene expression, matrix, Study,Breed,Sex,Age
}

}
})


write.csv(p_value_matrix,"p_statistic_114_tissues.matrix.csv")
write.csv(t_value_matrix,"t_statistic_114_tissues.matrix.csv")

stopImplicitCluster()




