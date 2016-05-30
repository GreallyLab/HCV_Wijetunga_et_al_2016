# Read in data with statistical modeling pvalues

meth_final<-read.table("meth_final_HCV_polyn_jstest.6.5.2014.bed")

# turn meth into bed file

meth_final[,5]<-meth_final[,1] 
meth_final[,1]<-meth_final[,2]
meth_final[,2]<-meth_final[,3]
meth_final[,3]<-meth_final[,3]+1

#remove NAs in statistical modeling p-value

meth_final<-meth_final[-which(is.na(meth_final[,20])),]

#  select significant models

meth_final<-meth_final[which(p.adjust(meth_final[,20], "fdr")<0.05),]

# select loci with significant differences of at least 15  between malignant and control 

meth_final<-meth_final[which(meth_final[,19]<0.05),]
meth_final<-meth_final[which(abs(meth_final[,9]-meth_final[,7])>.15),]

# separate the effects of increased and decreased methylation

up_meth<-meth_final[which(meth_final[,17]>0),]
down_meth<-meth_final[which(meth_final[,17]<0),]
up_meth<-up_meth[which(up_meth[,4]==11),]
down_meth<-down_meth[which(down_meth[,4]==11),] 

# TFBS enrichment 

# download all ENCODE HepG2 TFBS peak bed files into a directory
directory_with_TFBS_data<-""

vecup_TFBS<-sapply(dir(directory_with_TFBS_data,pattern="gz", full.names=T),function(i){
TFBS_intersection<-try(bedTools.2in(bed1=up_meth, bed2=read.table(i),opt.string="-u"), silent=T)
A<-ifelse(!is.null(ncol(TFBS_intersection)),
 {return(nrow(TFBS_intersection))},
 {1})
 return(A)}) 

vecdown_TFBS<-sapply(dir(directory_with_TFBS_data,pattern="gz", full.names=T),function(i){
TFBS_intersection<-try(bedTools.2in(bed1=down_meth, bed2=read.table(i),opt.string="-u"), silent=T)
A<-ifelse(!is.null(ncol(TFBS_intersection)),
 {return(nrow(TFBS_intersection))},
 {1})
return(A)}) 
 
names(vecup_TFBS)<-do.call(rbind,strsplit(do.call(rbind,strsplit(do.call(rbind,strsplit(names(vecup_TFBS),"/"))[,6],"\\."))[,1],"Hepg2"))[,2]
names(vecdown_TFBS)<-do.call(rbind,strsplit(do.call(rbind,strsplit(do.call(rbind,strsplit(names(vecdown_TFBS),"/"))[,6],"\\."))[,1],"Hepg2"))[,2]

vectotal_TFBS<-sapply(dir(directory_with_TFBS_data,pattern="gz", full.names=T),function(i){

 TFBS_intersection<-try(bedTools.2in(bed1=meth_final, bed2=read.table(i),opt.string="-u"), silent=T)
A<-ifelse(!is.null(ncol(TFBS_intersection)),
 {return(nrow(TFBS_intersection)+1)},
 {1})
 return(A)})

vecup_TFBS<-(vecup_TFBS/nrow(up_meth))/(vectotal_TFBS/nrow(meth_final))
vecup_TFBS[which(vecup_TFBS==0)]<-min(vecup_TFBS[-which(vecup_TFBS==0)])
vecdown_TFBS<-(vecdown_TFBS/nrow(down))/(vectotal_TFBS/nrow(meth_final))
vecdown TFBS[which(vecdown_TFBS==0)]<-min(vecdown_TFBS[-which(vecdown_TFBS==0)]) 

# make Figure 4

par(mar=c(10,3,3,3))

barplot(sort(log(vecup_TFBS),decreasing=T),las=2,ylim=c(-2.5, 1))
abline(h=0)
barplot(sort(log(vecdown_TFBS),decreasing=T),las=2,ylim=c(-2.5, 1)) 
abline(h=0)



# TFBS enrichment permutation tests

permute_TFBS_up<-sapply(directory_with_TFBS_data, function(i){
return(
bed_boot(
bed1=up_meth, 
bed2=read.table(i),
bed3=meth_final,
n.iter=100))
})

permute_TFBS_down<-sapply(directory_with_TFBS_data, function(i){
return(
bed_boot(
bed1=down_meth, 
bed2=read.table(i),
bed3=meth_final,
n.iter=100))
})

# plot TFBS permutation tests

par(mfrow=c(11,2),mar=c(1,2,1,1));
 for( i in seq(2,62,3)){plot(density(permute_TFBS_up[[i]]),xlim=c(0,120),main=i); abline(v=permute_TFBS_up[[i]])}

par(mfrow=c(11,2),mar=c(1,2,1,1));
 for( i in seq(2,62,3)){plot(density(permute_TFBS_down[[i]]),xlim=c(0,120),main=i); abline(v=permute_TFBS_down[[i]])}

