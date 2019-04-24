if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")

library(TCGAbiolinks)
query_crc_tmr<-GDCquery(project = 'TCGA-COAD',
                         legacy = TRUE,
                         data.category = 'Gene expression',
                         experimental.strategy = 'RNA-Seq',
                         platform = 'Illumina HiSeq',
                         data.type = "Gene expression quantification",
                        file.type = "results",
                sample.type = 'Primary solid Tumor'
                         )


GDCdownload(query_crc_tmr)
tcga_crc<-GDCprepare(query=query_crc_tmr,save=TRUE,save.filename='crc_tumor1')

View(tcga_crc)
query_crc_nor<-GDCquery(project = 'TCGA-COAD',
                        legacy = TRUE,
                        data.category = 'Gene expression',
                        experimental.strategy = 'RNA-Seq',
                        platform = 'Illumina HiSeq',
                        data.type = "Gene expression quantification",
                        file.type = "results",
                        sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor)
tcga_nor<-GDCprepare(query=query_crc_nor,save=TRUE,save.filename='crc_nor1')
View(tcga_nor)

query_crc_tmr2<-GDCquery(project = 'TCGA-READ',
                        legacy = TRUE,
                        data.category = 'Gene expression',
                        experimental.strategy = 'RNA-Seq',
                        platform = 'Illumina HiSeq',
                        data.type = "Gene expression quantification",
                        file.type = "results",
                        sample.type = 'Primary solid Tumor'
)

GDCdownload(query_crc_tmr2)
tcga_crc2<-GDCprepare(query=query_crc_tmr2,save=TRUE,save.filename='crc_tumor2')

View(tcga_crc)
query_crc_nor2<-GDCquery(project = 'TCGA-READ',
                        legacy = TRUE,
                        data.category = 'Gene expression',
                        experimental.strategy = 'RNA-Seq',
                        platform = 'Illumina HiSeq',
                        data.type = "Gene expression quantification",
                        file.type = "results",
                        sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor2)
tcga_nor2<-GDCprepare(query=query_crc_nor2,save=TRUE,save.filename='crc_nor2')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DelayedArray", version = "3.8")


library(SummarizedExperiment)
sum_coad<-SummarizedExperiment::cbind(tcga_crc,tcga_nor)
sum_read<-SummarizedExperiment::cbind(tcga_crc2,tcga_nor2)
sum_crc<-SummarizedExperiment::cbind(sum_coad,sum_read)


View(sum_crc)


data<- assay(sum_crc)
gene_info<-rowRanges(sum_crc)

View(sum_crc)

gene_dataf<-as.data.frame(data)
#prepare clinical data
clin_info<-sum_crc@colData
View(clin_info)

clin<-data.frame(clin_info@rownames,clin_info$patient,clin_info$gender,clin_info$year_of_birth,clin_info$year_of_death,clin_info$tumor_stage,clin_info$definition,clin_info$primary_diagnosis,clin_info$age_at_diagnosis)

sum_crc$vital_status


View(data)

library(dplyr) 
data<-as.data.frame(data)
data[2,2]



data=t(data)
gene_dataf=t(gene_dataf)
gene_dataf[1,0]


rownames(data[2,])
class(data)
write.table(gene_dataf,file='E:/final_thesis/data_from_TCGA/datao.csv',sep=',')
write.table(clin,file='E:/final_thesis/data_from_TCGA/clin.csv',sep=',')

#-------------------

# download the FPKM-uq values

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")

library(TCGAbiolinks)
query_crc_tmr<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - FPKM-UQ",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Primary solid Tumor'
)


GDCdownload(query_crc_tmr)
tcga_crc<-GDCprepare(query=query_crc_tmr,save=TRUE,save.filename='crc_tumor1')

View(tcga_crc)
query_crc_nor<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - FPKM-UQ",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor)
tcga_nor<-GDCprepare(query=query_crc_nor,save=TRUE,save.filename='crc_nor1')
View(tcga_nor)

query_crc_tmr2<-GDCquery(project = 'TCGA-READ',
                         data.category = 'Transcriptome Profiling',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = "HTSeq - FPKM-UQ",
                         data.type = "Gene Expression Quantification",
                         sample.type = 'Primary solid Tumor'
)

GDCdownload(query_crc_tmr2)
tcga_crc2<-GDCprepare(query=query_crc_tmr2,save=TRUE,save.filename='crc_tumor2')

View(tcga_crc)
query_crc_nor2<-GDCquery(project = 'TCGA-READ',
                         data.category = 'Transcriptome Profiling',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = "HTSeq - FPKM-UQ",
                         data.type = "Gene Expression Quantification",
                         sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor2)
tcga_nor2<-GDCprepare(query=query_crc_nor2,save=TRUE,save.filename='crc_nor2')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DelayedArray", version = "3.8")


library(SummarizedExperiment)
sum_coad<-SummarizedExperiment::cbind(tcga_crc,tcga_nor)
sum_read<-SummarizedExperiment::cbind(tcga_crc2,tcga_nor2)
sum_crc<-SummarizedExperiment::cbind(sum_coad,sum_read)


View(sum_crc)


data<- assay(sum_crc)
gene_info<-rowRanges(sum_crc)

View(sum_crc)

gene_dataf<-as.data.frame(data)
#prepare clinical data
clin_info<-sum_crc@colData
View(clin_info)

clin<-data.frame(clin_info@rownames,clin_info$patient,clin_info$gender,clin_info$year_of_birth,clin_info$year_of_death,clin_info$tumor_stage,clin_info$definition,clin_info$primary_diagnosis,clin_info$age_at_diagnosis)

View(data)

library(dplyr) 
data<-as.data.frame(data)
data[2,2]



data=t(data)
gene_dataf=t(gene_dataf)
gene_dataf[1,0]


rownames(data[2,])
class(data)
write.table(gene_dataf,file='E:/final_thesis/data_from_TCGA/data_n/data.csv',sep=',')
write.table(clin,file='E:/final_thesis/data_from_TCGA/data_n/clin.csv',sep=',')

#-------------
# download the FPKM values

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")

library(TCGAbiolinks)
query_crc_tmr<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - FPKM",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Primary solid Tumor'
)


GDCdownload(query_crc_tmr)
tcga_crc<-GDCprepare(query=query_crc_tmr,save=TRUE,save.filename='crc_tumor1')

#raw counts
query_crc_rowcounts<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - Counts",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Primary solid Tumor'
)


GDCdownload(query_crc_rowcounts)
tcga_crcrow<-GDCprepare(query=query_crc_rowcounts,save=TRUE,save.filename='crc_tumorrow')

gene_info<-rowRanges(tcga_crcrow)
data<- as.data.frame(data)

View(tcga_crc)
query_crc_nor<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - FPKM",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor)
tcga_nor<-GDCprepare(query=query_crc_nor,save=TRUE,save.filename='crc_nor1')
View(tcga_nor)

query_crc_tmr2<-GDCquery(project = 'TCGA-READ',
                         data.category = 'Transcriptome Profiling',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = "HTSeq - FPKM",
                         data.type = "Gene Expression Quantification",
                         sample.type = 'Primary solid Tumor'
)

GDCdownload(query_crc_tmr2)
tcga_crc2<-GDCprepare(query=query_crc_tmr2,save=TRUE,save.filename='crc_tumor2')

View(tcga_crc)
query_crc_nor2<-GDCquery(project = 'TCGA-READ',
                         data.category = 'Transcriptome Profiling',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = "HTSeq - FPKM",
                         data.type = "Gene Expression Quantification",
                         sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor2)
tcga_nor2<-GDCprepare(query=query_crc_nor2,save=TRUE,save.filename='crc_nor2')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DelayedArray", version = "3.8")


library(SummarizedExperiment)
sum_coad<-SummarizedExperiment::cbind(tcga_crc,tcga_nor)
sum_read<-SummarizedExperiment::cbind(tcga_crc2,tcga_nor2)
sum_crc<-SummarizedExperiment::cbind(sum_coad,sum_read)


View(sum_crc)


data<- assay(sum_crc)
gene_info<-rowRanges(sum_crc)

View(sum_crc)

gene_dataf<-as.data.frame(data)
#prepare clinical data
clin_info<-sum_crc@colData
View(clin_info)

clin<-data.frame(clin_info@rownames,clin_info$patient,clin_info$gender,clin_info$year_of_birth,clin_info$year_of_death,clin_info$tumor_stage,clin_info$definition,clin_info$primary_diagnosis,clin_info$age_at_diagnosis)

View(data)

library(dplyr) 
data<-as.data.frame(data)
data[2,2]



data=t(data)
gene_dataf=t(gene_dataf)
gene_dataf[1,0]


rownames(data[2,])
class(data)
write.table(gene_dataf,file='E:/final_thesis/data_from_TCGA/FPKMdata/data.csv',sep=',')
write.table(clin,file='E:/final_thesis/data_from_TCGA/FPKMdata/clin.csv',sep=',')


#______raw counts-----


if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")

library(TCGAbiolinks)
query_crc_tmr_rawcounts<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - Counts",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Primary solid Tumor'
)


GDCdownload(query_crc_tmr_rawcounts)
tcga_crc_counts<-GDCprepare(query=query_crc_tmr_rawcounts,save=TRUE,save.filename='crc_tumorc')


query_crc_nor_counts<-GDCquery(project = 'TCGA-COAD',
                        data.category = 'Transcriptome Profiling',
                        experimental.strategy = 'RNA-Seq',
                        workflow.type = "HTSeq - Counts",
                        data.type = "Gene Expression Quantification",
                        sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor_counts)
tcga_nor_counts<-GDCprepare(query=query_crc_nor_counts,save=TRUE,save.filename='crc_norc')
View(tcga_nor)

query_crc_tmr2_counts<-GDCquery(project = 'TCGA-READ',
                         data.category = 'Transcriptome Profiling',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = "HTSeq - Counts",
                         data.type = "Gene Expression Quantification",
                         sample.type = 'Primary solid Tumor'
)

GDCdownload(query_crc_tmr2_counts)
tcga_crc2_counts<-GDCprepare(query=query_crc_tmr2_counts,save=TRUE,save.filename='crc_tumor2c')


query_crc_nor2_count<-GDCquery(project = 'TCGA-READ',
                         data.category = 'Transcriptome Profiling',
                         experimental.strategy = 'RNA-Seq',
                         workflow.type = "HTSeq - Counts",
                         data.type = "Gene Expression Quantification",
                         sample.type = 'Solid Tissue Normal'
)

GDCdownload(query_crc_nor2_count)
tcga_nor2_counts<-GDCprepare(query=query_crc_nor2_count,save=TRUE,save.filename='crc_nor2c')


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DelayedArray", version = "3.8")


library(SummarizedExperiment)
sum_coad<-SummarizedExperiment::cbind(tcga_crc_counts,tcga_nor_counts)
sum_read<-SummarizedExperiment::cbind(tcga_crc2_counts,tcga_nor2_counts)
sum_crc<-SummarizedExperiment::cbind(sum_coad,sum_read)


View(sum_crc)


data<- assay(sum_crc)
gene_info<-rowRanges(sum_crc)

View(sum_crc)

gene_dataf<-as.data.frame(data)
#prepare clinical data
clin_info<-sum_crc@colData
View(clin_info)

clin<-data.frame(clin_info@rownames,clin_info$patient,clin_info$gender,clin_info$year_of_birth,clin_info$year_of_death,clin_info$tumor_stage,clin_info$definition,clin_info$primary_diagnosis,clin_info$age_at_diagnosis)

View(data)

library(dplyr) 
data<-as.data.frame(data)
data[2,2]

yh<-clin_info$vital_status=='dead'
length(yh[yh==TRUE])
clin_info$vital_status

clin_info$days_to_last_follow_up
clin_info$days_to_death



data=t(data)
gene_dataf=t(gene_dataf)
gene_dataf[1,0]


rownames(data[2,])
class(data)
write.table(gene_dataf,file='E:/final_thesis/data_from_TCGA/raw counts/data.csv',sep=',')
write.table(clin,file='E:/final_thesis/data_from_TCGA/raw counts/clin.csv',sep=',')
