# make figure S8

library("BSgenome.Hsapiens.UCSC.hg19")

cg_mat_h3k27me3<-uph3k4me1[uph3k4me1[,7]%in%Genes[Genes[,2]%in%cbind(setdiff(group2, group1)),1],][,1:3]
cg_mat_h3k27me3<-t(apply(cg_mat_h3k27me3, 1, function(i){
chr<-Hsapiens[[as.character(i[1])]]
active(masks(chr))<-F
windows<-Views(chr,  start=as.numeric(i[2]), end=as.numeric(i[3]))
a<-alphabetFrequency(windows,baseOnly=T)
cg<-dinucleotideFrequency(windows)[,7]
return(c(i, a ,cg))
}))
cg_mat_h3k27me3<-as.data.frame(cg_mat_h3k27me3)
cg_mat_h3k27me3<-cbind(cg_mat_h3k27me3,(as.numeric(cg_mat_h3k27me3[,9])/(as.numeric(cg_mat_h3k27me3[,5])*as.numeric(cg_mat_h3k27me3[,6])))*((as.numeric(cg_mat_h3k27me3[,3])-as.numeric(cg_mat_h3k27me3[,2]))-as.numeric(cg_mat_h3k27me3[,8])))

cg_mat_both<-uph3k4me1[uph3k4me1[,7]%in%Genes[Genes[,2]%in%cbind(intersect(group1, group2)),1],][,1:3]
cg_mat_both<-t(apply(cg_mat_both, 1, function(i){
chr<-Hsapiens[[as.character(i[1])]]
active(masks(chr))<-F
windows<-Views(chr,  start=as.numeric(i[2]), end=as.numeric(i[3]))
a<-alphabetFrequency(windows,baseOnly=T)
cg<-dinucleotideFrequency(windows)[,7]
return(c(i, a ,cg))
}))
cg_mat_both<-as.data.frame(cg_mat_both)
cg_mat_both<-cbind(cg_mat_both,(as.numeric(cg_mat_both[,9])/(as.numeric(cg_mat_both[,5])*as.numeric(cg_mat_both[,6])))*((as.numeric(cg_mat_both[,3])-as.numeric(cg_mat_both[,2]))-as.numeric(cg_mat_both[,8])))

cg_mat_h3k4me1<-uph3k4me1[uph3k4me1[,7]%in%Genes[Genes[,2]%in%cbind(setdiff(group1, group2)),1],][,1:3]
cg_mat_h3k4me1<-t(apply(cg_mat_h3k4me1, 1, function(i){
chr<-Hsapiens[[as.character(i[1])]]
active(masks(chr))<-F
windows<-Views(chr,  start=as.numeric(i[2]), end=as.numeric(i[3]))
a<-alphabetFrequency(windows,baseOnly=T)
cg<-dinucleotideFrequency(windows)[,7]
return(c(i, a ,cg))
}))
cg_mat_h3k4me1<-as.data.frame(cg_mat_h3k4me1)
cg_mat_h3k4me1<-cbind(cg_mat_h3k4me1,(as.numeric(cg_mat_h3k4me1[,9])/(as.numeric(cg_mat_h3k4me1[,5])*as.numeric(cg_mat_h3k4me1[,6])))*((as.numeric(cg_mat_h3k4me1[,3])-as.numeric(cg_mat_h3k4me1[,2]))-as.numeric(cg_mat_h3k4me1[,8])))

 plot.multi.dens <- function(s)
{
junk.x = NULL
junk.y = NULL
for(i in 1:length(s))
{
junk.x = c(junk.x, density(s[[i]])$x)
junk.y = c(junk.y, density(s[[i]])$y)
}
xr <- range(junk.x)
yr <- range(junk.y)
plot(density(s[[1]]), xlim = xr, ylim = yr, main = "",)
for(i in 1:length(s))
{
lines(density(s[[i]]), xlim = xr, ylim = yr, col = i, lwd=3)

}
}

plot.multi.dens(list(cg_mat_h3k4me1[,10],cg_mat_h3k27me3[,10],cg_mat_both[,10]))

