#for ORIGINAL permutation
fe_2<-read.csv('cancer.csv')
fe_1<-read.csv('normal.csv')


fe_2$gender<-NULL

fe_1$gender<-NULL

#install.packages("corrplot")
library(corrplot)


corr<-cor(fe_2, method = "spearman")
corrplot(corr, method="color", type = "lower",tl.cex = 0.6,number.cex = 0.5,
         tl.col="black",tl.srt=45, sig.level = 0.05,insig = "blank")


corr<-cor(fe_1, method = "spearman")
corrplot(corr, method="color", type = "lower",tl.cex = 0.6,number.cex = 0.5,
         tl.col="black",tl.srt=45, sig.level = 0.05,insig = "blank")





#----------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", version = "3.8")

#library(edgeR)
#library(limma)


k<-read.csv('22dec.csv')

#prepare the experimental design matrix
#k$clin_info.gender<-NULL

k$outcome<-as.factor(k$outcome)
levels(k$outcome)

#group<-factor(k$outcome, levels = c('0','1'))
#design<-model.matrix(~group)
#colnames(design)<-c('Cancer','Normal vs Cancer')

#design

#prepare the so called EList object.
counts<-k[,-1]
counts<-t(counts)
counts
#dge<-DGEList(counts=counts, group = group)
#keep <- filterByExpr(dge, design)
#dge <- dge[keep,,keep.lib.sizes=FALSE]
dim(k)

#v <- voom(dge,design, plot=TRUE)

#fit <- lmFit(v, design)
#fit <- eBayes(fit)
#topTable(fit, coef='Normal vs Cancer',number = Inf,sort.by="M")


#topTable(fit, coef=ncol(design),number = Inf,sort.by="P")


#------ use DEseq 
library(DESeq2)
library(apeglm)

cold<-factor(k$outcome)
dds<-DESeqDataSetFromMatrix(counts, DataFrame(cold),~cold)


nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

dds <- DESeq(dds)


# for calculate the rlog___


vst=varianceStabilizingTransformation(dds)  
vsd=assay(vst)  # vsd is now the normalized log2-transformed data
save(dds,vsd,cold,file="rlogm.RData") # save it for the future
vsd

rlg=t(vsd)
rlg=as.data.frame(rlg)
rlg=cbind(cold,rlg)
rlg


write.table(rlg,file = 'rlogmmmix.csv',sep=',')

# hierarchical clustering of samples and heatmap of sample similarities
library(pheatmap)
pheatmap(cor(vsd))  # yes it is that simple

# Principle coordiante analysis
library(vegan)
library(rgl)
library(ape)

colData(dds)$cold=='0'


colData(dds)$'0'
# assembling table of conditions to lable PCoA plot:
# (in the chunk below, replace factor1 and factor2 with your actual factor names from myConditions table)
factor1=as.character(colData(dds)$cold=='0')
factor2=as.character(colData(dds)$cold=='1')
oneByTwo=paste(factor1,factor2,sep=".")
conditions=data.frame(cbind(factor1,factor2,oneByTwo))

# actual PCoA analysis
dds.pcoa=pcoa(vegdist(t(vsd),method="manhattan")/1000)
scores=dds.pcoa$vectors

# plotting
plot(scores[,1], scores[,2],col=as.numeric(as.factor(factor1)))
ordispider(scores,factor2,label=T)
ordiellipse(scores,factor2)

# interactive 3d plot - can rotate it by dragging mouse
radiusScale=2 # change it if the spheres come out too small or too big in the next one
plot3d(scores[,1], scores[,2],scores[,3],col=as.numeric(as.factor(factor1)),type="s",radius=radiusScale*as.numeric(as.factor(factor2)))

# formal permutation-based analysis of variance 
adonis(t(vsd)~factor1*factor2,data=conditions,method="manhattan")  

#______

res <- results(dds)

res 
showres <- results(dds,tidy = TRUE)

summary(res)
resultsNames(dds)

resLFC <- lfcShrink(dds, coef=2, type="apeglm")

summary(resLFC)

male_log<-as.data.frame(res)


write.table(female_log,file='male_log.csv', sep = ',')


library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

res = res[order(res$pvalue),]
#row.names(res)

showres
showres$symbol = mapIds(org.Hs.eg.db,
                    keys=showres$row, 
                    column="SYMBOL",
                    keytype="ENSEMBL",
                    multiVals="first")
showres$entrez = mapIds(org.Hs.eg.db,
                    keys=showres$row, 
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")
showres$name =   mapIds(org.Hs.eg.db,
                    keys=showres$row, 
                    column="GENENAME",
                    keytype="ENSEMBL",
                    multiVals="first")


list1<-res$symbol
list1
res

write.table(list1, file = 'list52.csv')





#___________________



k<-read.csv('47counts.csv')
k$outcome
k$outcome<-as.factor(k$outcome)
levels(k$outcome)
counts<-counts[,-1]
dim(counts)

scaled_counts<-scale(counts)


write.table(scaled_counts, file="malenormalized_counts.csv", sep=",")
#_____


counts<-t(counts)
counts
dim(k)


library(DESeq2)
MedianNorm(counts, alternative=TRUE)



outcome<-factor(k$outcome)
dds<-DESeqDataSetFromMatrix(counts, DataFrame(outcome),~outcome)



dds <- estimateSizeFactors( dds )
sizeFactors(dds)


logcounts <- log2( counts(dds, normalized=TRUE) + 1 )
logcounts

f<-normalizeMedianValues(logcounts)
f

normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts<-t(normalized_counts)
normalized_counts
write.table(normalized_counts, file="fenormalized_counts.txt", sep="\t", quote=F)
#______________________________________________________________________
#


library(purrr)
library(dplyr)

k<-read.csv('female.csv')

k$ENSG00000000971
#for male
k<-k[!(k$outcome == 0),]
k$outcome
k$outcome<-NULL
k$clin_info.gender<-NULL

#for female
k$clin_info.definition

k<-k[!(k$clin_info.definition == 0),]

k$clin_info.definition<-NULL

k$clin_info.gender<-NULL
k$patient_age<-NULL


counts<-k


dim(counts)

for (i in 1: ncol(counts)) {
  counts[[i]]<-as.numeric(counts[[i]])
}


scaled_counts<-scale(counts,center = TRUE)



scaled_counts<-as.data.frame(scaled_counts)
scaled_counts$ENSG00000183134
dim(scaled_counts)


####_____
fe_s<-read.csv('female_s.csv')

fe_sur<-cbind(fe_s,scaled_counts)
dim(fe_sur)


fe_sur$status

typeof(fe_sur$ENSG00000122641)


fe_sur<-fe_sur[!is.na(fe_sur$status),]
fe_sur<-fe_sur[!is.na(fe_sur$survival_days),]



write.table(fe_sur, file="female_sum_sur.csv", sep=",")






library(dplyr)
dplyr::filter(fe_sur,  !is.na(columnname))







write.table(f, file="male_sum_sur.csv", sep=",")


m_s<-read.csv('male_s.csv')
m<-read.csv('manormalized_counts.csv')
m_s$status
m$ENSG00000000003
m_sur<-cbind(m_s,m)

write.table(m_sur, file="male_sum_sur.csv", sep=",")




