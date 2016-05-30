# make figure S5
bcelldata<-read.table("/Users/Ari/Desktop/GreallyLab/HCV/Final.BCell.polyn.6.26.14.txt")
bcelldata2<-bcelldata[-which(is.na(bcelldata[,14])),]
plot(bcelldata2[,12], -log(bcelldata2[,13]),xlim=c(-2,2),pch=19,cex=.5,ylim=c(0,5),col=alpha("blue",.5))