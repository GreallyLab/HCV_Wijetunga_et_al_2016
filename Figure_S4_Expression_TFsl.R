# Make figure S4

c_index<-which(do.call(rbind,strsplit(colnames(expression)[2:28],"_"))[,2]=="C")
nt_index<-which(do.call(rbind,strsplit(colnames(expression)[2:28],"_"))[,2]=="NT")
t_index<-which(do.call(rbind,strsplit(colnames(expression)[2:28],"_"))[,2]=="T")
par(mfrow=c(2,4))
#Expression of FOXA2
vioplot(log(as.numeric(expression[grep("NM_021784", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_021784", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_021784", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))
#Expression of FoxA1
vioplot(log(as.numeric(expression[grep("NM_004496", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_004496", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_004496", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))
#Expression of HNF4A
vioplot(log(as.numeric(expression[grep("NM_001258355", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_001258355", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_001258355", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))
#Expression of MAFF
vioplot(log(as.numeric(expression[grep("NM_001161572", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_001161572", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_001161572", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))
#Expression of MAFK
vioplot(log(as.numeric(expression[grep("NM_002360", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_002360", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_002360", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))
 #Expression of CEBPb
vioplot(log(as.numeric(expression[grep("NM_005194", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_005194", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_005194", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))
#Expression of RXRA
vioplot(log(as.numeric(expression[grep("NM_002957", expression[,1]),][2:28][c_index])+1),log(as.numeric(expression[grep("NM_002957", expression[,1]),][2:28][nt_index])+1),log(as.numeric(expression[grep("NM_002957", expression[,1]),][2:28][t_index])+1), col="gray",border=NA, rectCol="gray29", wex=.5,colMed="red",pchMed=16, ylim=c(0,7))

