
# for integrating clinical and transcriptomic information, and divide into female and male corhorts. 
library(dplyr)
clin<-read.csv('clin.csv')
gene_data<-read.csv('data.csv')



sum_data<-cbind(clin,gene_data)
sum_data$id <- NULL
write.table(sum_data,file='E:/final_thesis/data_from_TCGA/raw counts/sum_data.csv',sep=',')


sum_data$clin_info.age_at_diagnosis <- as.integer(sum$clin_info.age_at_diagnosis/365)
sum_data$clin_info.age_at_diagnosis

colnames(sum_data)

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

