# crctranscriptomeproject_ml

This project is Biomarker Discovery and Identification for Colorectal Cancer Regarding Gender Differences with Machine Learning Techniques,
which is a master thesis written by YanghongLiu 

Organization:  Karolinska Institutet, Stockholm University, Scilifelab

In this study, the gene expression variances regarding gender discrepancy in colorectal cancer is explored by a series analysis 
including correlation analysis, principle component analysis, machine learning modeling, cancer hallmarks analysis and survival analysis. 

Codes for TCGA CRC transcriptome download, feature selection, machine learning, correlation analysis, PCA and survival analysis will be provided.

Tools used in this study: R, Python

note:
1. data preprocessing details including merging clinic and transcriptome information, seperating female and male cohort, converting rlog value, median normalization for survival analysis an so on will not be provided. Please email me if you have any question: genevieveyanghong@gmail.com
2. PCA will be included in machine learningxxx.ipynb
3. the size of the data downloaded from TCGA is too large, the dataset for female and male cohorts are compressed and stored in the 'fpkm' file under 'data' file
