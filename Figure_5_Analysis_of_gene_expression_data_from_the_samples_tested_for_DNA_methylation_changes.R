
# Read in expression data 
expression<-read.table("Merged_Gene_Counts_31_samples.DESeq.normalized.pvalues.batchcorrected.txt",header=T)

#select only significantly differentially expressed genes

expression_sig<-expression_sig[which(expression_sig$polyn_lrt_p<0.005),]
expression_sig<-expression_sig[which( expression_sig$polyn_p_TvsC<0.05),]

# select upregulated and downregulated genes

upexp<-expression_sig[which(expression_sig$polyn_B_TvsC>0),]
downexp<-expression_sig[which(expression_sig$polyn_B_TvsC<0),]

# convert gene names from refseq to gene symbols

genes<-read.table("Gene_to_Refseq.txt")
upexpgenes<-unique(genes[genes[,1]%in%upexp[,1],2])
downexpgenes<-unique(genes[genes[,1]%in%downexp[,1],2])


# make figure S1 top

kmeansdata<-read.table("Means_Post_Combat_for_Kmeans_10.13.14.txt", header=T)

# select only means

kmeansdata_means<-(kmeansdata[,7:9])

# Divide each mean by the control group and take log transformation

for(i in 3:1){ kmeansdata_means[,i]<-kmeansdata_means[,i]/kmeansdata_means[,1]}
kmeansdata_means<-log(kmeansdata_means)

# find WSS and plot figure s1 bottom

wss <- (nrow(kmeansdata_means)-1)*sum(apply(kmeansdata_means,2,var))
for (i in 1:20) { print(i); wss[i] <- sum(kmeans(kmeansdata_means,
                                      centers=i)$withinss)}
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

# make figure  S1 bottom
 
kmeansdata<-read.table("Means_Post_Combat_for_Kmeans_10.13.14.exp.txt", header=F)
kmeansdata<-exp(kmeansdata[,-1])

# Divide each mean by the control group and take log transformation

for(i in 3:1){ kmeansdata[,i]<-kmeansdata[,i]/kmeansdata[,1]}
kmeansdata_means<-log(kmeansdata)

# find WSS

wss <- (nrow(mydata)-1)*sum(apply(kmeansdata_means,2,var))
  for (i in 1:20) { print(i); wss[i] <- sum(kmeans(kmeansdata_means,
                                       centers=i)$withinss)}
plot(1:20, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# make figure 5a     
# find kmeans clustering with 6 centers

kmeans_exp<-kmeans(kmeansdata, centers=6)
 
# find confidence intervals for clusters with bootstrapping
 
inputdata<-kmeansdata_means
clustering<-kmeans_exp

error<-sapply(1:6, function(i){return(sapply(2:3,function(j){print(i); print(j);
in.data<-inputdata[which(clustering$cluster==i),j];
dist_out<-sapply(1:100,function(x){return(mean(in.data[order(rnorm(length(in.data)))[1:100]]))})
return(c(quantile(dist_out,.025),quantile(dist_out,.975)))
}))})
error<-t(error)
error<-rbind(cbind(0, error[,1], error[,3]),cbind(0, error[,2], error[,4]))
error<-error[c(1,7,2,8,3,9,4,10,5,11,6,12),]

par(mfrow=c(3,2))
for(k in 1:6){
 plot.new()
plot.window(xlim=c(1,3),ylim=c(-1.5,1.5))
temp<-inputdata[which(clustering$cluster==k),]
apply(temp[order(rnorm(nrow(temp)))[1:100],], 1, function(x){lines(1:3, x,col="gray",lwd=.5)})
polygon(c(1:3,3:1),c(error[k*2-1,], rev(error[k*2,])), col="blue", border=NA)
 lines(1:3, clustering$centers[k,] ,col="red")
abline(h=0)
text(2, 0,table(clustering$cluster)[k])
axis(2)  
axis(1)
}

# make figure 5b

plot(log2(expression$C_mean+1), log2(expression$T_mean+1), xlim=c(0,20),ylim=c(0,20),pch=16,col="gray")
points(log2(upexp$C_mean+1), log2(upexp$T_mean+1), col="red")
points(log2(downexp$C_mean+1), log2(downexp$T_mean+1), col="green")

# Make figure 5c

library(venneuler)
plot(venneuler(c(A=622,B=105,C=291,"B&C"=0,"B&A"=4,"C&A"=18,"B&C&A"=0)))


# make figure 5d

library(beanplot)

expression_TCGA<-read.table("TCGA_TNT_EXP_HCV.txt",header=T)
tcgamethorig<-read.table("TCGA_TNT_Meth_HCV.bed",header=F)


nt_expdown<-expression_TCGA[expression_TCGA[,1]%in%downexpgenes,3]
t_expdown<-expression_TCGA[expression_TCGA[,1]%in%downexpgenes,4]
nt_expup<-expression_TCGA[expression_TCGA[,1]%in%upexpgenes,3]
t_expup<-expression_TCGA[expression_TCGA[,1]%in%upexpgenes,4]

par(mfrow=c(2,1))
beanplot(as.numeric(nt_expup), as.numeric(t_expup),col="red",wex=.75,boxwex=.75,what=c(0,1,1,1),method="overplot",beanlines="median",ll=.05,beanlinewd=.5,border=F,ylim=c(0,15.2))
beanplot(as.numeric(nt_expdown), as.numeric(t_expdown),col="lightgreen",wex=.75,boxwex=.75,what=c(0,1,1,1),method="overplot",beanlines="median",ll=.05,beanlinewd=.5,border=F,ylim=c(0,15.2))

uph3k4me1<-
bedTools.2in(bed1=up_meth,bed2=
bedTools.2in(bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"), 
bed1=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),functionstring="subtractBed"),opt.string="-wa -wb")
uph3k4me1<-uph3k4me1[,23:25]
uph3k4me1<-bedTools.2in(functionstring="closestBed", bed1=uph3k4me1, bed2=read.table("hg19_refseqgenes_symbol_tss.bed"), opt.string="-d")
uph3k4me1<-uph3k4me1[which(abs(uph3k4me1[,10])<5000),]
uph3k4me1tcga<-bedTools.2in(bed1=tcgamethorig,bed2=uph3k4me1,opt.string="-wa -wb")

uph3k4me1tcga<-uph3k4me1tcga[-which(duplicated(uph3k4me1tcga[,4])),]
beanplot(uph3k4me1tcga[,6], uph3k4me1tcga[,7],col="lightblue",wex=.75,boxwex=.75,what=c(0,1,1,1),method="overplot",beanlines="median",ll=.05,beanlinewd=.5,border=F)
