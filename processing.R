
setwd('D:/R_data/tcgabio')

clin<-read.csv('clin.csv')
View(clin)
gene_data<-read.csv('data.csv')

library(dplyr)

sum_data<-cbind(clin,gene_data)
sum_data$id <- NULL
write.table(sum_data,file='E:/final_thesis/data_from_TCGA/raw counts/sum_data.csv',sep=',')
sum_data$


sum_data$clin_info.age_at_diagnosis <- as.integer(sum$clin_info.age_at_diagnosis/365)
sum_data$clin_info.age_at_diagnosis

colnames(sum_data)
sum_data<-cbind(clin,gene_data)

sum_data$clin_info.year_of_death<-NULL
sum_data$clin_info.rownames<- NULL
sum_data$clin_info.patient <-NULL
sum_data$clin_info.year_of_birth<-NULL
sum_data$clin_info.tumor_stage <- NULL
sum_data$clin_info.primary_diagnosis <- NULL

sum_data$clin_info.age_at_diagnosis<-NULL

sum_data$clin_info.definition<-as.character(sum_data$clin_info.definition)

for (i in 1:nrow(sum_data)){
  if (sum_data$clin_info.definition[i]=='Primary solid Tumor')
    sum_data$clin_info.definition[i]<-'1'
  else if(sum_data$clin_info.definition[i]=='Solid Tissue Normal')
    sum_data$clin_info.definition[i]<-'0'
}
sum_data$clin_info.definition<-as.factor(sum_data$clin_info.definition)



colnames(sum_data)[colnames(sum_data)=="clin_info.definition"] <- "outcome"



female_sum <- subset(sum_data, clin_info.gender=='female')
male_sum <- subset(sum_data, clin_info.gender=='male')

write.table(female_sum,file='E:/final_thesis/data_from_TCGA/raw counts/female.csv',sep=',')
write.table(male_sum,file='E:/final_thesis/data_from_TCGA/raw counts/male.csv',sep=',')

#______ counts match with the 80

mc<-read.csv('final.csv')
mc$ENSG00000000003
mc$num<-NULL
mc$ID<-NULL
mc$Gender<-NULL
mc$Status<-as.factor(mc$Status)

mc=t(mc)
mc<-as.data.frame(mc)

write.table(mc,file='E:/final_thesis/data_from_TCGA/raw counts/mc.csv',sep=',')

mc$ENSG00000000003
male_sum$clin_info.gender<-NULL


# read in counts file 
mcounts<-read.csv('mc.csv')


for (i in 2: ncol(fe_data)) {
  fe_data[[i]]<-as.numeric(fe_data[[i]])
}

## match the list with counts values
mc$num<-NULL
mc$ID<-NULL
mc$Gender<-NULL

colnames(mc)[colnames(mc)=='Status']<-'outcome'

female_sum$ENSG00000000003

fff<-read.csv('47DE.csv')
class(fff)
fff<-t(fff)

for (i in 1: ncol(fff)) {
  fff[[i]]<-as.character(fff[[i]])
}


fff<-as.vector(fff)
fff



selected_fe1<- subset(male_sum,select=c(fff))

selected_fe1<-cbind(male_sum$outcome,selected_fe1)


names(selected_fe1)[names(selected_fe1) == "male_sum$outcome"] <- "outcome"


write.table(selected_fe1,file='47counts.csv',sep=',')


