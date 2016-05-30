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


# Find background distribution

 vector_meth_final<-c(

 nrow(bedTools.2in(bed1=bedTools.2in(
bed1=meth_final, 
bed2=read.table("hg19_cgislands.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgshores.bed"),opt.string="-u"))
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=meth_final, 
bed2=read.table("hg19_cgshores.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgislands.bed"),opt.string="-u")))
,

length(unique(bedTools.2in(bed1=bedTools.2in(
bed1=meth_final, 
bed2=read.table("hg19_cgshores.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgislands.bed"),opt.string="-v")[,2]))
 ,
 
 nrow((bedTools.2in(
 bed1=meth_final, 
 bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-u")))  
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=meth_final, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-v"), 
bed2=read.table("hg19_refseqgenes.bed"),opt.string="-u")))
,

length(unique(bedTools.2in(bed1=bedTools.2in(
bed1=meth_final, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-v"), 
bed2=read.table("hg19_refseqgenes.bed"),opt.string="-v")[,2]))
,

nrow((bedTools.2in(
bed1=meth_final,
bed2=read.table("H3K4me3.Adult_Liver.signal.meth_final.bed"),opt.string="-u")))
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=meth_final,
bed2=read.table("H3K4me3.Adult_Liver.signal.meth_final.bed"),opt.string="-v"), 
bed2=read.table("H3K4me1.Adult_Liver.signal.meth_final.bed"),opt.string="-u")))
,

nrow((bedTools.2in(
bed1=meth_final,
bed2=read.table("H3K36me3.Adult_Liver.signal.meth_final.bed"),opt.string="-u")))
,

nrow((bedTools.2in(bed1=meth_final,
bed2=read.table("H3K27me3.Adult_Liver.signal.meth_final.bed"),opt.string="-u")))
,
nrow(bedTools.2in(bed1=(bedTools.2in(bed1=bedTools.2in(
bed1=meth_final,
bed2=read.table("H3K4me3.Adult_Liver.signal.meth_final.bed"),opt.string="-v"), 
bed2=read.table("H3K4me1.Adult_Liver.signal.meth_final.bed"),opt.string="-u")),
bed2=read.table("H3K27ac.Adult_Liver.signal.meth_final.bed"),opt.string="-u"))


) 

# Enrichment for increasing and decreasing loci

vecup<-c(
 
nrow((
bedTools.2in(bed1=bedTools.2in(
bed1=up_meth, 
bed2=read.table("hg19_cgislands.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgshores.bed"),opt.string="-u")))
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=up_meth, 
bed2=read.table("hg19_cgshores.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgislands.bed"),opt.string="-u")))
,

length(unique(bedTools.2in(bed1=bedTools.2in(
bed1=up_meth, 
bed2=read.table("hg19_cgshores.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgislands.bed"),opt.string="-v")[,2]))
,

nrow((bedTools.2in(
bed1=up_meth, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-u")))  
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=up_meth, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-v"), 
bed2=read.table("hg19_refseqgenes.bed"),opt.string="-u")))
,

length(unique(bedTools.2in(bed1=bedTools.2in(
bed1=up_meth, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-v"), 
bed2=read.table("hg19_refseqgenes.bed"),opt.string="-v")[,2]))
,

nrow((bedTools.2in(
bed1=up_meth,
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=up_meth,
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"),opt.string="-v"), 
bed2=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow((bedTools.2in(
bed1=up_meth,
bed2=read.table("H3K36me3.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow((bedTools.2in(
bed1=up_meth,
bed2=read.table("H3K27me3.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow(bedTools.2in(bed1=(bedTools.2in(bed1=bedTools.2in(
bed1=up_meth,
bed2=read.table("H3K4me3.Adult_Liver.signal.meth_final.bed"),opt.string="-v"), 
bed2=read.table("H3K4me1.Adult_Liver.signal.meth_final.bed"),opt.string="-u")),
bed2=read.table("H3K27ac.Adult_Liver.signal.meth_final.bed"),opt.string="-u"))

)



vecdown<-c(
 
nrow((bedTools.2in(bed1=bedTools.2in(
bed1=down_meth, 
bed2=read.table("hg19_cgislands.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgshores.bed"),opt.string="-u")))
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=down_meth, 
bed2=read.table("hg19_cgshores.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgislands.bed"),opt.string="-u")))
,

length(unique(bedTools.2in(bed1=bedTools.2in(
bed1=down_meth, 
bed2=read.table("hg19_cgshores.bed"),opt.string="-v"), 
bed2=read.table("hg19_cgislands.bed"),opt.string="-v")[,2]))
,

nrow((bedTools.2in(
bed1=down_meth, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-u")))  
,

nrow((bedTools.2in(bed1=bedTools.2in(
bed1=down_meth, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-v"), 
bed2=read.table("hg19_refseqgenes.bed"),opt.string="-u")))
,

length(unique(bedTools.2in(bed1=bedTools.2in(
bed1=down_meth, 
bed2=read.table("hg19_refseqgenes_tss_4kb.bed"), opt.string="-v"), 
bed2=read.table("hg19_refseqgenes.bed"),opt.string="-v")[,2]))
,

nrow((bedTools.2in(
bed1=down_meth,
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow((
bedTools.2in(bed1=bedTools.2in(
bed1=down_meth,
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"),opt.string="-v"), 
bed2=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow((bedTools.2in(
bed1=down_meth,
bed2=read.table("H3K36me3.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,

nrow((bedTools.2in(
bed1=down_meth,
bed2=read.table("H3K27me3.Adult_Liver.signal.FINAL.bed"),opt.string="-u")))
,
nrow(bedTools.2in(bed1=(bedTools.2in(bed1=bedTools.2in(
bed1=down_meth,
bed2=read.table("H3K4me3.Adult_Liver.signal.meth_final.bed"),opt.string="-v"), 
bed2=read.table("H3K4me1.Adult_Liver.signal.meth_final.bed"),opt.string="-u")),
bed2=read.table("H3K27ac.Adult_Liver.signal.meth_final.bed"),opt.string="-u"))

)

up_enrichment<-(vecup/nrow(up_meth))/( vector_meth_final/nrow(meth_final)) 
names(up_enrichment)<-c("CpG_Shores","CpG_Islands","CpG_Oceans","Ref_Promoters","Ref_Bodies","Intergenic","H3K4me3","H3K4me1(w/o me3)","H3K36me3","H3K27me3", "H3K27ac/H3K4me3 (w/o me1)")

down_enrichment<-(vecdown/nrow(down_meth))/( vector_meth_final/nrow(meth_final)) 
names(down_enrichment)<-c("CpG_Shores","CpG_Islands","CpG_Oceans","Ref_Promoters","Ref_Bodies","Intergenic","H3K4me3","H3K4me1(w/o me3)","H3K36me3","HK27me3","H3K27ac/H3K4me3 (w/o me1)")

# make Figure 3 

par(mfrow=c(2,1))
par(mar=c(5.1,4.1,1,2.1))
barplot(sort(log(up_enrichment),decreasing=T),ylim=c(-2.25,1.5),las=2,col=1,axes=F,axisnames=F,ylab="ln(O/E)")
abline(h=0)
axis(2)
text(names(sort(up_enrichment,decreasing=T)),x=seq(.7, 90.7, length.out=76), y=rep(-1.5
,76),srt=90,col=2)
barplot(sort(log(down_enrichment),decreasing=T),ylim=c(-2.25,1.5),las=2,col=1,axes=F,axisnames=F,ylab="ln(O/E)")
axis(2)
text(names(sort(down_enrichment,decreasing=T)),x=seq(.7, 90.7, length.out=76), y=rep(-1.5
,76),srt=90,col=2)
abline(h=0)

# permutation tests

promoters<-read.table("hg19_refseqgenes_symbol_tss_4KB.bed")

 bodies<-bedTools.2in(bed1=read.table("hg19_refseqgenes_symbol.bed"),
bed2=read.table("hg19_refseqgenes_symbol_tss_4KB.bed"), functionstring="subtractBed")

intergenic<-bedTools.2in(bed1=bedTools.2in(bed1=read.table("hg19.bed"),bed2=read.table("hg19_refseqgenes_symbol.bed"),
functionstring="subtractBed"), bed2=read.table("hg19_refseqgenes_symbol_tss_4KB.bed"),
functionstring="subtractBed")

 cgislands<-bedTools.2in(bed1=read.table("hg19_cgislands.bed"),
bed2=read.table("hg19_cgshores.bed"), functionstring="subtractBed")

cgshores<-bedTools.2in(bed2=read.table("hg19_cgislands.bed"),
bed1=read.table("hg19_cgshores.bed"), functionstring="subtractBed")

cgoceans<-bedTools.2in(bed1=bedTools.2in(bed1=read.table("hg19.bed"),bed2=read.table("hg19_cgislands.bed"),
functionstring="subtractBed"), bed2=read.table("hg19_cgshores.bed"),
functionstring="subtractBed")


H3K4me1<-bedTools.2in(bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"), 
bed1=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),functionstring="subtractBed")

H3K36me3<-read.table("H3K36me3.Adult_Liver.signal.FINAL.bed")

H3K27me3<-read.table("H3K27me3.Adult_Liver.signal.FINAL.bed")

H3K4me3<-read.table("H3K4me3.Adult_Liver.signal.FINAL.bed")



promoters_annot_up<-bed_boot(bed1=up_meth, bed2="hg19_refseqgenes_tss_4kb.bed", bed3=meth_final)

bodies_annot_up<-bed_boot(bed1=up_meth, bed2=
bedTools.2in(bed1=read.table("hg19_refseqgenes_symbol.bed"),
bed2=read.table("hg19_refseqgenes_symbol_tss_4KB.bed"), functionstring="subtractBed")
, bed3=meth_final)

intergenic_annot_up<-bed_boot(bed1=up_meth, bed2=bedTools.2in(bed1=bedTools.2in(bed1=read.table("hg19.bed"),bed2=read.table("hg19_refseqgenes_symbol.bed"),
functionstring="subtractBed"), bed2=read.table("hg19_refseqgenes_symbol_tss_4KB.bed"),
functionstring="subtractBed"), bed3=meth_final)

shores_annot_up<-bed_boot(bed1=up_meth, bed2=bedTools.2in(bed2=read.table("hg19_cgislands.bed"),
bed1=read.table("hg19_cgshores.bed"), functionstring="subtractBed"), bed3=meth_final)

islands_annot_up<-bed_boot(bed1=up_meth, bed2=bedTools.2in(bed1=read.table("hg19_cgislands.bed"),
bed2=read.table("hg19_cgshores.bed"), functionstring="subtractBed"), bed3=meth_final)

oceans_annot_up<-bed_boot(bed1=up_meth, bed2=bedTools.2in(bed1=bedTools.2in(bed1=read.table("hg19.bed"),bed2=read.table("hg19_cgislands.bed"),
functionstring="subtractBed"), bed2=read.table("hg19_cgshores.bed"),
functionstring="subtractBed")
, bed3=meth_final)

h3k4me1_annot_up<-bed_boot(bed1=up_meth, bed2=bedTools.2in(bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"), 
bed1=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),functionstring="subtractBed")
, bed3=meth_final)

h3k36me3_annot_up<-bed_boot(bed1=up_meth, bed2=read.table("H3K36me3.Adult_Liver.signal.FINAL.bed"), bed3=meth_final)

h3k27me3_annot_up<-bed_boot(bed1=up_meth, bed2=read.table("H3K27me3.Adult_Liver.signal.FINAL.bed"), bed3=meth_final)

h3k4me3_annot_up<-bed_boot(bed1=up_meth, bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed")
, bed3=meth_final)

h3k27ac_annot_up<-bed_boot(bed1=up_meth, 
bed2=bedTools.2in(bed1=bedTools.2in(bed1=read.table("H3K27ac.Adult_Liver.signal.FINAL.bed"), 
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed")), 
bed2=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),functionstring="subtractBed")
, bed3=meth_final)




promoters_annot_down<-bed_boot(bed1=down_meth, bed2=promoters, bed3=meth_final)

bodies_annot_down<-bed_boot(bed1=down_meth, bed2=bodies, bed3=meth_final)

intergenic_annot_down<-bed_boot(bed1=down_meth, bed2=intergenic, bed3=meth_final)

shores_annot_down<-bed_boot(bed1=down_meth, bed2=cgshores, bed3=meth_final)

islands_annot_down<-bed_boot(bed1=down_meth, bed2=cgislands, bed3=meth_final)

oceans_annot_down<-bed_boot(bed1=down_meth, bed2=cgoceans, bed3=meth_final)

h3k4me1_annot_down<-bed_boot(bed1=down_meth, bed2=H3K4me1, bed3=meth_final)

h3k36me3_annot_down<-bed_boot(bed1=down_meth, bed2=H3K36me3, bed3=meth_final)

h3k27me3_annot_down<-bed_boot(bed1=down_meth, bed2=H3K27me3, bed3=meth_final)

h3k4me3_annot_down<-bed_boot(bed1=down_meth, bed2=H3K4me3, bed3=meth_final)

h3k27ac_annot_down<-bed_boot(bed1=down_meth, 
bed2=bedTools.2in(bed1=bedTools.2in(bed1=read.table("H3K27ac.Adult_Liver.signal.FINAL.bed"), 
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed")), 
bed2=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),functionstring="subtractBed")
, bed3=meth_final)

# plot enrichment permutation tests
par(mfrow=c(10,1)) 
plot(density(log(h3k4me1_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k4me1_annot_up[[1]]))
plot(density(log(h3k36me3_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k36me3_annot_up[[1]]))
plot(density(log(bodies_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(bodies_annot_up[[1]]))
plot(density(log(shores_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(shores_annot_up[[1]]))
plot(density(log(oceans_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(oceans_annot_up[[1]]))
plot(density(log(h3k4me3_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k4me3_annot_up[[1]]))
plot(density(log(intergenic_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(intergenic_annot_up[[1]]))
plot(density(log(promoters_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(promoters_annot_up[[1]]))
plot(density(log(islands_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(islands_annot_up[[1]]))
plot(density(log(h3k27me3_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k27me3_annot_up[[1]]))
plot(density(log(h3k27ac_annot_up[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k27ac_annot_up[[1]]))
 
par(mfrow=c(10,1)) 
plot(density(log(h3k27me3_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k27me3_annot_down[[1]]))
plot(density(log(intergenic_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(intergenic_annot_down[[1]]))
plot(density(log(oceans_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(oceans_annot_down[[1]]))
plot(density(log(shores_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(shores_annot_down[[1]]))
plot(density(log(bodies_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(bodies_annot_down[[1]]))
plot(density(log(islands_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(islands_annot_down[[1]]))
plot(density(log(promoters_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(promoters_annot_down[[1]]))
plot(density(log(h3k4me3_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k4me3_annot_down[[1]]))
plot(density(log(h3k36me3_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k36me3_annot_down[[1]]))
plot(density(log(h3k4me1_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k4me1_annot_down[[1]]))
plot(density(log(h3k27ac_annot_down[[2]])), xlim=c(1.5,7),main="",xlab="",ylab="")
abline(v=log(h3k27ac_annot_down[[1]]))


