# this code was running for three different cohorts, data processing part is slightly different, user can adjust based on needs
library(dplyr)
#the input is the three cohorts of female, male and mixed-sex FPKM value csv file
fe_data <- read.csv('female.csv')

#the processing might be different, it depends on the csv file you are dealing with , the colname will be different
fe_data$clin_info.gender<-NULL
fe_data$patient_age<-NULL

names(fe_data)[names(fe_data) == "clin_info.definition"] <- "outcome"

fe$outcome<-as.factor(fe$outcome)


for (i in 2: ncol(fe_data)) {
  fe_data[[i]]<-as.numeric(fe_data[[i]])
}


#feature selection with Vita

library(ranger)
library(vita)
library("randomForest") 

cv_vi = CVPVI(fe_data[,-1],fe_data[,1],k = 2,mtry = 252,ntree = 1000, ncores = 2)
print(cv_vi)

#compare them with the original permutation variable importance

cl.rf = randomForest(fe_data[,-1],fe_data[,1], mtry = 252,ntree = 500, importance = TRUE)
print(cl.rf)
round(cbind(importance(cl.rf, type=1, scale=FALSE),cv_vi$cv_varim),digits=15)


cv_p = NTA(cv_vi$cv_varim)
summary(cv_p,pless = 0.001)

#take the results from CVPVI, and store it into a file, then we have the ranking list with P-value 
fe=cv_p$pvalue

write.table(fe,file='vita_final.csv',sep=',')


## do this after vita, manually pick up those gene with P-values = 0 in excel file, and store them in .csv file(in txt file and conver it into csv format)

fff<-read.csv('list.csv')
class(fff)
fff<-t(fff)

for (i in 1: ncol(fff)) {
  fff[[i]]<-as.character(fff[[i]])
}


fff<-as.vector(fff)

# map the selected gene with their FPKM value

selected_fe1<- subset(fe,select=c(fff))


selected_fe1<-cbind(fe$outcome,selected_fe1)

names(selected_fe1)[names(selected_fe1) == "fe$outcome"] <- "outcome"

write.table(selected_fe1,file='33fpkm.csv',sep=',')


dim(selected_fe1)


# for selecting 1000 attributes using Boruta
library(ranger)
library(Boruta)

aa<-read.csv('33fpkm.csv')


aa$outcome

aa$outcome<-as.factor(aa$outcome)

for (i in 2: ncol(aa)) {
  aa[[i]]<-as.numeric(aa[[i]])
}


set.seed(1234)
boruta.aa<-Boruta(outcome~.,data = aa, doTrace = 2)
print(boruta.aa)

set.seed(1234)

selected_fe1$outcome<-as.factor(selected_fe1$outcome)
boruta.fe<-Boruta(outcome~.,data=selected_fe1, doTrace = 2)
print(boruta.fe)

#show genes confirmed as important, then we got the list
for (i in 1:length(boruta.fe$finalDecision)){
  if (boruta.fe$finalDecision[i]=='Confirmed')
    print(boruta.fe$finalDecision[i])
}


## do this after boruta, keep the list and then map the gene names with their fpkm values
fff<-read.csv('boruta71.csv')
class(fff)
fff<-t(fff)

for (i in 1: ncol(fff)) {
  fff[[i]]<-as.character(fff[[i]])
}


fff<-as.vector(fff)
fff

selected_fe2<- subset(fe_data,select=c(fff))

selected_fe2<-cbind(fe_data$outcome,selected_fe2)

names(selected_fe2)[names(selected_fe2) == "fe_data$outcome"] <- "outcome"

write.table(selected_fe2,file='fe_boruta2.csv',sep=',')




#  use mRMR after Vita
library(devtools)
#install_github("bhklab/mRMRe")
library(mRMRe)

selected_fe1$outcome<-ordered(selected_fe1$outcome)


f<- mRMR.data(data = data.frame(target=selected_fe1$outcome, selected_fe1[,-1]))
results <- mRMR.classic("mRMRe.Filter", data=f,target_indices = 1, feature_count = 100)


a <- as.data.frame(solutions(results))
print(solutions(results))

mrmr_100<-print(apply(solutions(results)[[1]], 2, function(x, y) { return(y[x]) }, y=featureNames(f)))
mrmr_100<-as.vector(mrmr_100)


se_mrmr_100<- subset(fe_data,select=c(mrmr_100))


se_mrmr_100<-cbind(fe_data$outcome,se_mrmr_100)

names(se_mrmr_100)[names(se_mrmr_100) == "fe_data$outcome"] <- "outcome"

write.table(se_mrmr_100,file='male_mrmr.csv',sep=',')


