# crctranscriptomeproject_ml

This project is Biomarker Discovery and Identification for Colorectal Cancer Regarding Sex Differences with Machine Learning Techniques by Yanghong Liu.

Organization:  Scilifelab, Stockholm University, Karolinska Institute

In this study, the gene expression variances regarding sex discrepancy in colorectal cancer is explored by a series analysis 
including correlation analysis, principle component analysis(PCA), machine learning modelling and survival analysis. 

Codes for TCGA CRC transcriptome download, feature selection, machine learning modelling, correlation analysis, PCA and survival analysis will be provided.

Tools used in this study: R, Python

datasets: (please notice FPKM and raw counts of the transcriptome were used in this study)
TCGA-female cohort, TCGA-male cohort, Swedish mixed-gender cohort 


Feature selection: There are three different algorithms used for feature selection: Vita, Boruta and mRMR. They were combined to optimize the feature selection results. 


NOTE:
1. Data preprocessing detail for merging clinic and transcriptome information and dividing into female and male cohorts are provided in data_processing.R under 'Data retrieval and pre-processing file'. Since the TCGA dataset is too large, even the compressed file is not able to be uploaded in github. The processed female and male datasets are not provided here. 
2. Other data processing details including DEGs analysis with DESeq2, converting into rlog value for ML, median normalization for survival analysis will not be provided. Please email me if you have any question on this: genevieveyanghong@gmail.com.
3. PCA, resampling and importance ranking will be included in machine learningxxx.ipynb
4.Processed datasets for machine learning and survival analysis are in file 'rlog for ML' and 'normalized dataset for survival analysis' under data file.
5. Please adjust the code based on your own case. You might have some changes based on which kind of transcriptome value you are using. This project is only for guidance.
