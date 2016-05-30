# read in adjusted methylation data

methdata<-read.table("Combat_batch_adjusted_logit.HCV.txt",header=T)

# read in annotation

HCV_annotation<-read.table("HCV_Tissue_HCV_annotationation.txt",header=T,sep="\t")

# select only methylation values

onlymethdata<-methdata[,5:ncol(methdata)]
onlymethdata<-onlymethdata[, order(colnames(onlymethdata))]

# make figure 1

# find most variable loci
methylation_variation<-apply(onlymethdata, 1, var)
for(i in 1:ncol(onlymethdata)){onlymethdata[,i]<-exp(onlymethdata[,i])/
(1+exp(onlymethdata[,i]))
}

TissueCols<-HCV_annotation$Tissue
TissueCols[which(TissueCols=="T")]<-"black"
TissueCols[which(TissueCols=="NT")]<-"gray"
TissueCols[which(TissueCols=="C")]<-"white"
onlymethdata<-onlymethdata[order(methylation_variation,decreasing=T),]

library(GMD)

heatmap.3(onlymethdata[1:5000,],Rowv=T,Colv=T,dendrogram="both",scale="none",breaks=seq(0,1, length.out=101),ColSideColors=as.matrix(TissueCols),col=colorRampPalette(c("yellow","blue"))(100),vclust.method="complete")
