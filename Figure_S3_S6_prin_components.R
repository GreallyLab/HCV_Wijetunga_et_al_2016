methdata<-read.table("HELPdata.cauchy.MM.totals.10.13.13.txt",header=T)

expdata<-read.table("Merged_Gene_Counts_31_samples.DESeq.normalized.txt",header=T)

Annot<-read.table("HCV_Tissue_Annotation.txt",header=T,sep="\t")

onlymethdata<-methdata[,5:ncol(methdata)]
onlymethdata<-onlymethdata[, order(colnames(onlymethdata))]

Annot_meth<-as.data.frame(Annot[Annot[,1]%in%colnames(onlymethdata),])
Annot_meth<-Annot_meth[,-19]
Annot_meth<-Annot_meth[,c(1,2,12,5,21,6,7,8,9,10,11,20,13,14,15,16,17,18,4,3,19)]


expdata<-expdata[,order(colnames(expdata))]
Annot_exp<-as.data.frame(Annot[Annot[,1]%in%colnames(expdata),])
Annot_exp<-Annot_exp[,-18]
Annot_exp<-Annot_exp[,c(1,2,12,5,21,6,7,8,9,10,11,20,13,14,15,16,17,18,4,3,19)]

library(stats)
princomp_meth=princomp(onlymethdata)
index<-c("factor","factor","continuous","continuous","continuous","continuous","continuous","continuous", "continuous","continuous","continuous","continuous", "continuous","continuous","continuous","continuous", "factor","continuous","factor","factor")
pvalmat_meth<-t(sapply(2:ncol(Annot_meth),function(i){
print(i)
if(index[i-1]=="factor"){return(sapply(1:10,function(j){
return(anova(lm(princomp_meth$loadings[,j]~as.factor(Annot_meth[,i])))[[5]][1])}))}

if(index[i-1]=="continuous"){return(sapply(1:10,function(j){
return(anova(lm(princomp_meth$loadings[,j]~(as.numeric(Annot_meth[,i]))))[[5]][1])}))}

}
))


princomp_exp=princomp(expdata)
index<-c("factor","factor","continuous","continuous","continuous","continuous","continuous","continuous", "continuous","continuous","continuous","continuous", "continuous","continuous","continuous","continuous", "factor","continuous","factor","factor")

pvalmat_exp<-t(sapply(2:ncol(Annot_exp),function(i){
print(i)
if(index[i-1]=="factor"){return(sapply(1:10,function(j){
return(anova(lm(princomp_exp$loadings[,j]~as.factor(Annot_exp[,i])))[[5]][1])}))}

if(index[i-1]=="continuous"){return(sapply(1:10,function(j){
return(anova(lm(princomp_exp$loadings[,j]~as.numeric(Annot_exp[,i])))[[5]][1])}))}

}
))

