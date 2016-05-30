# permutation test function

bed_boot<-function(bed1, bed2,bed3, n.iter=100){

mainobs<-try(bedTools.2in(bed1=bed1, bed2=bed2), silent=T)
A<-ifelse(!is.null(ncol(mainobs)),
 {(length(which(!duplicated(paste(mainobs[,1],mainobs[,2])))))},
 {1})
mainobs<-A
 
simulations<-sapply(1:n.iter, function(x){
print(x)
newbed<-bed3[order(rnorm(1:nrow(bed3)))[1:nrow(bed1)],]
temp<-try(bedTools.2in(bed1=newbed, bed2=bed2), silent=T)
ifelse(!is.null(ncol(temp)),
 {return(length(which(!duplicated(temp[,2]))))},
 {1})})

wilcox.p<-wilcox.test(simulations, mu=mainobs)$p.value

return(list(mainobs, simulations, wilcox.p))}
