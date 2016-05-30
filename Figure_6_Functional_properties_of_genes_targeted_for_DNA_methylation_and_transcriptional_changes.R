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

genes<-read.table("Gene_to_Refseq.txt")

uph3k4me1<-
bedTools.2in(bed1=up_meth,bed2=
bedTools.2in(bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"), 
bed1=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),functionstring="subtractBed"),opt.string="-wa -wb")
uph3k4me1<-uph3k4me1[,23:25]
uph3k4me1<-bedTools.2in(functionstring="closestBed", bed1=uph3k4me1, bed2=read.table("hg19_refseqgenes_symbol_tss.bed"), opt.string="-d")
uph3k4me1<-uph3k4me1[which(abs(uph3k4me1[,10])<5000),]

expression<-read.table("Merged_Gene_Counts_31_samples.DESeq.normalized.pvalues.batchcorrected.txt",header=T)


# make figure 6a

# H3K4me1 in liver to h3k4me1 in hepg2

group11<-cbind(unique(Genes[Genes[,1]%in%
bedTools.2in(bed1=bedTools.2in(
bed1=uph3k4me1, 
bed2=read.table("H3K27me3.Adult_Liver.signal.FINAL.bed"), opt.string="-v"), 
bed2=read.table("wgEncodeBroadHistoneHepg2H3k04me1StdAlnRep0.bam_VS_wgEncodeBroadHistoneHepg2ControlStdAlnRep0.bam_peaks.gappedPeak.gz"), opt.string="-wa -wb")[,7],2]))

# H3K4me1 in liver to H3K27me3 in hepg2

group2<-cbind(unique(Genes[Genes[,1]%in%
bedTools.2in(bed1=bedTools.2in(
bed1=uph3k4me1, 
bed2=read.table("H3K27me3.Adult_Liver.signal.FINAL.bed"), opt.string="-v"), 
bed2=read.table("wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep0.bam_VS_wgEncodeBroadHistoneHepg2ControlStdAlnRep0.bam_peaks.broadPeak.gz"), opt.string="-wa -wb")[,7],2]))

# genes with both h3k4me1 and h3k27me3 in hepg2

cbind(intersect(group1, group2))

# genes with one or the other

cbind(setdiff(group1, group2))
cbind(setdiff(group2, group1))

# make figure 6b

options(stringsAsFactors=F)

#to install

library(SMITE)

#load methylation p values
methylation<-meth_final

methylation<-methylation[,c(1:3, 17,19)]
methylation<-replace_zeros(methylation, 5)
methylation<-methylation[-which(is.na(methylation[,5])),]
 
#load gene expression p values
genes<-expression

genes<-cbind(rownames(genes), genes)
genes<-genes[,c(1,33,29, 36)]

genes[,1]<-convertGeneIds(gene_IDs=genes[,1], ID_type ="refseq", ID_convert_to ="symbol")
genes<-genes[-which(is.na(genes[,1])),]

genes<-split(genes, genes[,1])
genes<-lapply(genes, function(i){if(nrow(as.data.frame(i))>1)
{i<-i[which(i[,4]==min(i[,4],na.rm=T))[1],]}; return(i)}) 
genes<-do.call(rbind, genes)
genes<-genes[,-1]
genes<-replace(genes[, 3], genes[,3] == 0, min(subset(genes[, 3], genes[, 3] != 0), na.rm=TRUE))
expression<-genes
colnames(expression)<-c("effect","jstat","pval")

#Create annotation with gene symbols and liver h3k4me1

data(hg19_genes_bed)

h3k4me1<-bedTools.2in(
bed2=read.table("H3K4me3.Adult_Liver.signal.FINAL.bed"), 
bed1=read.table("H3K4me1.Adult_Liver.signal.FINAL.bed"),
functionstring="subtractBed")

HCV_annotation<-makePvalueAnnotation( data=hg19_genes, otherdata=list(h3k4me1=h3k4me1), gene_name_col=5,other_tss_distance=5000)

#fill in expression data
HCV_annotation<-annotateExpression(HCV_annotation, expression, effect_col=1, pval_col=3)

#fill in methylation data
#this step takes ~10 minutes
HCV_annotation<-annotateModification(HCV_annotation, methylation, weight_by_method=c(promoter="distance",body="distance",h3k4me1="distance"),verbose=T,mod_corr=T)

#create a pvalue object that will count the effect of the h3k4me1 as bidirectional
HCV_annotation<-makePvalueObject(HCV_annotation, effect_directions=c(methylation_promoter="decrease",methylation_body="increase",methylation_h3k4me1="bidirectional"))

#normalize the pvalues compared to colExp
HCV_annotation<-normalizePval(HCV_annotation,ref="expression",method="rescale")
	
#score with all four features contributing
HCV_annotation<-scorePval(HCV_annotation,weights=c(methylation_promoter=.3,methylation_body=.1,expression=.6,methylation_h3k4me1=.6))

#load REACTOME 
load(system.file("data","Reactome.Symbol.Igraph.rda", package="SMITE"))
 
#run epimods using REACTOME network
HCV_annotation<-runSpinglass(HCV_annotation, REACTOME,maxsize=50, num_interations=1000, simplify=T)

plotModule (HCV_annotation, which.network=1) 
