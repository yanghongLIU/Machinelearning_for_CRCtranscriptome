#for ORIGINAL permutation
fe_1<-read.csv('cancer.csv')
fe_2<-read.csv('normal.csv')


fe_1$gender<-NULL

fe_2$gender<-NULL

library(corrplot)


corr<-cor(fe_1, method = "spearman")
corrplot(corr, method="color", type = "lower",tl.cex = 0.6,number.cex = 0.5,
         tl.col="black",tl.srt=45, sig.level = 0.05,insig = "blank")


corr<-cor(fe_2, method = "spearman")
corrplot(corr, method="color", type = "lower",tl.cex = 0.6,number.cex = 0.5,
         tl.col="black",tl.srt=45, sig.level = 0.05,insig = "blank")







