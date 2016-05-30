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


plot(density(c(up_meth[,17], down_meth[,17])))


# Read in means data

kmeansdata<-read.table("Means_Post_Combat_for_Kmeans_10.13.14.txt", header=T)

# select only means

kmeansdata_means<-(kmeansdata[,7:9])

# Divide each mean by the control group and take log transformation

for(i in 3:1){ kmeansdata_means[,i]<-kmeansdata_means[,i]/kmeansdata_means[,1]}
kmeansdata_means<-log(kmeansdata_means)

# find WSS and plot figure s1 top

wss <- (nrow(kmeansdata_means)-1)*sum(apply(kmeansdata_means,2,var))
for (i in 1:20) { print(i); wss[i] <- sum(kmeans(kmeansdata_means,
                                      centers=i)$withinss)}
plot(1:20, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")


# make figure 2b

# find kmeans clustering with 8 centers

kmeans_meth<-kmeans(kmeansdata_means, centers=8)

# find confidence intervals for clusters with bootstrapping

inputdata<-kmeansdata_means
clustering<-kmeans_meth
error<-sapply(1:8, function(i){return(sapply(2:3,function(j){print(i); print(j);
in.data<-inputdata[which(clustering$cluster==i),j];
dist_out<-sapply(1:n,function(x){return(mean(in.data[order(rnorm(length(in.data)))[1:100]]))})
return(c(quantile(dist_out,.025),quantile(dist_out,.975)))
}))})
error<-t(error)
error<-rbind(cbind(0, error[,1], error[,3]),cbind(0, error[,2], error[,4]))
error<-error[c(1,9,2,10,3,11,4,12,5,13,6,14,7,15,8,16),]
 
 
par(mfrow=c(3,3))
for(k in 1:8){
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
barplot(table(kmeans_meth$cluster)[order(kmeans_meth$centers[,3])],horiz=T)
text(rep(400000, 8), 1:8, table(kmeans_meth$cluster)[order(kmeans_meth$centers[,3])])
