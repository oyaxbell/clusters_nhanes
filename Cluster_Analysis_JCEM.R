#Trends in the prevalence of diabetes subgroups in U.S adults: A data-driven cluster analysis in NHANES from 1988 to 2018 
#Analysis: Neftali Eduardo Antonio-Villa (nefoantonio@hotmail.com); Luisa Fernandez-Chirino (fernandez.luisa@comunidad.unam.mx); Omar Yaxmehen Bello-Chavolla (oyaxbell@yahoo.com.mx)
#Diclosure: All data is available from: https://www.cdc.gov/nchs/nhanes/index.htm

#----Library loading----
library(haven); 
library(tidyverse); library(dplyr); 
library(survival)
library(ggthemes)
library(ggsci)
library(survey)
library(segmented)
library(ggpubr)
library(quantreg)
library(data.table)
library(nhanesA)

library(tidyverse) #for all data wrangling
library(cowplot) #for manuscript ready figures
library(lme4) #for lmer & glmer models
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats)
library(ggplot2)

#----Dataset managment-----
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/PROYECTOS/T2D CLUSTERS")
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/T2D CLUSTERS")

NHANES_DM2<-readr::read_rds("NHANES_DM2_III_2019_2.rds")


#Include diabetes 1999-2004 (2000-2004)

DIQ.2000<-nhanes("DIQ")%>%dplyr::select(SEQN,DIQ040Q)%>%rename(seqn=SEQN,Age_diabetes=DIQ040Q)
DIQ.2002<-nhanes("DIQ_B")%>%dplyr::select(SEQN,DID040Q)%>%rename(seqn=SEQN,Age_diabetes=DID040Q)
DIQ.2004<-nhanes("DIQ_C")%>%dplyr::select(SEQN,DID040Q)%>%rename(seqn=SEQN,Age_diabetes=DID040Q)
DIQ_2000_2004<-rbind(DIQ.2000, DIQ.2002, DIQ.2004)
DIQ_2000_2004$seqn<-paste0(DIQ_2000_2004$seqn,"-4")

NHANES_DM2<-left_join(NHANES_DM2, DIQ_2000_2004, by = c("seqn")) %>%
  mutate(Age_diabetes.y=as.numeric(Age_diabetes.y),
         Age_diabetes.x=as.numeric(Age_diabetes.x))
NHANES_DM2$Age_diabetes.y[NHANES_DM2$Age_diabetes.y==99999]<-NA  
NHANES_DM2<-NHANES_DM2%>%mutate(Age_diabetes = coalesce(Age_diabetes.x, Age_diabetes.y)) %>% dplyr::select(-c(Age_diabetes.x,Age_diabetes.y))


#Descriptive Analysis
nrow(NHANES_DM2)
table(NHANES_DM2$AGE>=18)

#Dataset managment
NHANES_DM2$Glucose[NHANES_DM2$Glucose == 88888] <- NA
NHANES_DM2$Glycohemoglobin[NHANES_DM2$Glycohemoglobin == 8888] <- NA
NHANES_DM2$Insulin[NHANES_DM2$Insulin == 888888] <- NA
NHANES_DM2$time_diabetes<-(NHANES_DM2$AGE-NHANES_DM2$Age_diabetes)
table(is.na(NHANES_DM2$time_diabetes),NHANES_DM2$YEAR)

NHANES_DM2$DM2_final<-as.numeric((NHANES_DM2$Glycohemoglobin>=6.5 |NHANES_DM2$Glucose>=126 | NHANES_DM2$Dr_diabetes==1))
NHANES_DM2$DM2_final<-na.tools::na.replace(NHANES_DM2$DM2_final,0)

table(NHANES_DM2$DM2_final,NHANES_DM2$AGE>=18)
sum(!is.na(NHANES_DM2[NHANES_DM2$AGE_18==1,]$Cluster))
table(NHANES_DM2$YEAR,NHANES_DM2$Cluster)
prop.table(table(NHANES_DM2$YEAR,NHANES_DM2$Cluster),1)*100

#Variable recodification

NHANES_DM2$AGE_18<-NULL;NHANES_DM2$AGE_18[NHANES_DM2$AGE>=18]<-1;NHANES_DM2$AGE_18[NHANES_DM2$AGE<18]<-0
NHANES_DM2$BMI[NHANES_DM2$BMI == 8888] <- NA
NHANES_DM2$Weight[NHANES_DM2$Weight == 888888] <- NA
NHANES_DM2$Height[NHANES_DM2$Height == 88888.0] <- NA
NHANES_DM2$BMI_2<-NHANES_DM2$Weight/((NHANES_DM2$Height/100)^2)
NHANES_DM2<-NHANES_DM2%>%mutate(one = 1)

d1<-dummies::dummy(NHANES_DM2$Cluster)
NHANES_DM2<-cbind(NHANES_DM2,d1)


NHANES_DM2$YEAR_2<-NULL;
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="1988-1990"]<-1
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="1991-1994"]<-2
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="1999-2000"]<-3
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2001-2002"]<-3
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2003-2004"]<-4
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2005-2006"]<-4
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2007-2008"]<-5
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2009-2010"]<-5
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2011-2012"]<-6
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2013-2014"]<-6
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2015-2016"]<-7
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2017-2018"]<-7
NHANES_DM2$YEAR_2<-factor(NHANES_DM2$YEAR_2, labels = c("1988-1990","1991-1994","1999-2002","2003-2006","2007-2010","2011-2014","2015-2018"))



NHANES_DM2$RACE_REC<-NULL;
NHANES_DM2$RACE_REC[NHANES_DM2$RACE==1]<-1
NHANES_DM2$RACE_REC[NHANES_DM2$RACE==2]<-1
NHANES_DM2$RACE_REC[NHANES_DM2$RACE==3]<-2
NHANES_DM2$RACE_REC[NHANES_DM2$RACE==4]<-3
NHANES_DM2$RACE_REC[NHANES_DM2$RACE==5]<-2
NHANES_DM2$RACE_REC<-factor(NHANES_DM2$RACE_REC, labels = c("Mexican American","Non-Hispanic White", "Non-Hispanic Black"))


NHANES_DM2$age_cat<-NULL;
NHANES_DM2$age_cat[NHANES_DM2$AGE>=0 & NHANES_DM2$AGE<=17]<-0
NHANES_DM2$age_cat[NHANES_DM2$AGE>=18 & NHANES_DM2$AGE<=39]<-1
NHANES_DM2$age_cat[NHANES_DM2$AGE>=40 & NHANES_DM2$AGE<=64]<-2
NHANES_DM2$age_cat[NHANES_DM2$AGE>=65]<-3
NHANES_DM2$age_cat<-factor(NHANES_DM2$age_cat,labels = c("<18 years", "≥18-39 years",
                                                         "40-64 years","≥65 years"))

NHANES_DM2$OBESITY_2<-NULL;
NHANES_DM2$OBESITY_2[NHANES_DM2$BMI_2<30]<-0
NHANES_DM2$OBESITY_2[NHANES_DM2$BMI_2>=30]<-1
NHANES_DM2$OBESITY_2<-factor(NHANES_DM2$OBESITY_2, labels = c("Non-Obese", "Obese"))


NHANES_DM2$BMI_CAT<-NULL;
NHANES_DM2$BMI_CAT[NHANES_DM2$BMI_2<25]<-0
NHANES_DM2$BMI_CAT[NHANES_DM2$BMI_2>=25 & NHANES_DM2$BMI_2<30]<-1
NHANES_DM2$BMI_CAT[NHANES_DM2$BMI_2>=30]<-2
NHANES_DM2$BMI_CAT<-factor(NHANES_DM2$BMI_CAT, labels = c("Normal-Weight", "Overweight", "Obese"))


table(NHANES_DM2$YEAR,NHANES_DM2$BMI_CAT,NHANES_DM2$DM2_final,useNA = "always")

NHANES_DM2$GENDER_REC<-NULL;
NHANES_DM2$GENDER_REC[NHANES_DM2$GENDER==1]<-1
NHANES_DM2$GENDER_REC[NHANES_DM2$GENDER==2]<-0
NHANES_DM2$GENDER_REC<-factor(NHANES_DM2$GENDER_REC, labels = c("Women", "Men"))

#Education dataset

nhanes3<-read_csv("Clean Datasets/nhanes3.csv")
nhanes3_edu<-nhanes3%>%dplyr::select(SEQN,SDPPHASE,HFA8R)
nhanes3_edu$seqn<-paste0(str_pad(nhanes3_edu$SEQN, 1,pad = "0"),c("-3"))

nhanes3_edu$EDU_REC<-NULL
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R<=8]<-1
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R>=9 & nhanes3_edu$HFA8R<=11]<-2
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R==12]<-3
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R>=13 & nhanes3_edu$HFA8R<=16]<-4
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R==17]<-5
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R==88]<-9
nhanes3_edu$EDU_REC[nhanes3_edu$HFA8R==99]<-1
nhanes3_edu_2<-nhanes3_edu%>%dplyr::select(c(seqn,EDU_REC))

edu_n<-read_csv("Clean Datasets/educ.csv")
edu_n$seqn<-paste0(str_pad(edu_n$SEQN, 1,pad = "0"),c("-4"))
edu_n$EDU_REC<-edu_n$Schooling

edu_n_2<-edu_n%>%dplyr::select(seqn,EDU_REC)
edu_final<-rbind(nhanes3_edu_2,edu_n_2)

NHANES_DM2<-NHANES_DM2%>%left_join(edu_final,by="seqn")
NHANES_DM2$EDU_REC_2<-NULL
NHANES_DM2$EDU_REC_2[is.na(NHANES_DM2$EDU_REC)==1]<-1
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==1]<-1
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==2]<-2
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==3]<-2
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==4]<-3
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==5]<-3
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==7]<-1
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==9]<-1

NHANES_DM2$EDU_REC_2<-factor(NHANES_DM2$EDU_REC_2, labels = c("Primary School \nor \nNo-Education", "Secondary, High-School \nor \nAA Degree", "College \nor \nHigher"))


#Diabetes Recently diagnosed

NHANES_DM2$time_diabetes<-NHANES_DM2$AGE-NHANES_DM2$Age_diabetes
NHANES_DM2$time_diabetes[NHANES_DM2$time_diabetes>=0 & NHANES_DM2$DM2_final==0]<-NA
NHANES_DM2$time_diabetes[is.na(NHANES_DM2$time_diabetes) & NHANES_DM2$DM2_final==1]<-0

NHANES_DM2$RECENT_DIABETES<-NULL;
NHANES_DM2$RECENT_DIABETES[is.na(NHANES_DM2$time_diabetes)]<-1
NHANES_DM2$RECENT_DIABETES[NHANES_DM2$time_diabetes<5]<-2
NHANES_DM2$RECENT_DIABETES[NHANES_DM2$time_diabetes>=5]<-3
NHANES_DM2$RECENT_DIABETES<-as.factor(NHANES_DM2$RECENT_DIABETES)

#Recent MARD

NHANES_DM2$RECENT_MARD<-NULL;
NHANES_DM2$RECENT_MARD[NHANES_DM2$RECENT_DIABETES==1 & NHANES_DM2$ClusterMARD==0]<-1
NHANES_DM2$RECENT_MARD[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterMARD==0]<-2
NHANES_DM2$RECENT_MARD[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterMARD==0]<-3
NHANES_DM2$RECENT_MARD[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterMARD==1]<-4
NHANES_DM2$RECENT_MARD[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterMARD==1]<-5
NHANES_DM2$RECENT_MARD<-as.factor(NHANES_DM2$RECENT_MARD)

#Recent MOD

NHANES_DM2$RECENT_MOD<-NULL;
NHANES_DM2$RECENT_MOD[NHANES_DM2$RECENT_DIABETES==1 & NHANES_DM2$ClusterMOD==0]<-1
NHANES_DM2$RECENT_MOD[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterMOD==0]<-2
NHANES_DM2$RECENT_MOD[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterMOD==0]<-3
NHANES_DM2$RECENT_MOD[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterMOD==1]<-4
NHANES_DM2$RECENT_MOD[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterMOD==1]<-5
NHANES_DM2$RECENT_MOD<-as.factor(NHANES_DM2$RECENT_MOD)

#Recent SIRD

NHANES_DM2$RECENT_SIRD<-NULL;
NHANES_DM2$RECENT_SIRD[NHANES_DM2$RECENT_DIABETES==1 & NHANES_DM2$ClusterSIRD==0]<-1
NHANES_DM2$RECENT_SIRD[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterSIRD==0]<-2
NHANES_DM2$RECENT_SIRD[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterSIRD==0]<-3
NHANES_DM2$RECENT_SIRD[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterSIRD==1]<-4
NHANES_DM2$RECENT_SIRD[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterSIRD==1]<-5
NHANES_DM2$RECENT_SIRD<-as.factor(NHANES_DM2$RECENT_SIRD)

#Recent SIID

NHANES_DM2$RECENT_SIID<-NULL;
NHANES_DM2$RECENT_SIID[NHANES_DM2$RECENT_DIABETES==1 & NHANES_DM2$ClusterSIID==0]<-1
NHANES_DM2$RECENT_SIID[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterSIID==0]<-2
NHANES_DM2$RECENT_SIID[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterSIID==0]<-3
NHANES_DM2$RECENT_SIID[NHANES_DM2$RECENT_DIABETES==2 & NHANES_DM2$ClusterSIID==1]<-4
NHANES_DM2$RECENT_SIID[NHANES_DM2$RECENT_DIABETES==3 & NHANES_DM2$ClusterSIID==1]<-5
NHANES_DM2$RECENT_SIID<-as.factor(NHANES_DM2$RECENT_SIID)

#----Diabetes weighted subset----

NHANES_all <- svydesign(data=NHANES_DM2, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
NHANES_subset <- subset(NHANES_all, AGE_18==1)
NHANES_subset1<- subset(NHANES_all, AGE_18==1 & !is.na(NHANES_DM2$Insulin))

#----All Population characteristics (Table 1)-----

NHANES_DM2_Diab<-NHANES_DM2%>%dplyr::filter(AGE_18==1)%>%dplyr::filter(DM2_final==1)%>%dplyr::filter(!is.na(Cluster))
columns <- c('Parameter',"n = 5,489")

## Sex Female

sexo <- table(NHANES_DM2_Diab$GENDER_REC)[c(1)]
sexoprop <- round(prop.table(table(NHANES_DM2_Diab$GENDER_REC,useNA = "always")),4)[c(1)]*100

Sexo<-`names<-`(as.data.frame(matrix(c("Female (%)",paste(sexo,paste0('(',sexoprop,')'))),ncol=2,byrow=T)),columns)

#Age


Edad<-`names<-`(as.data.frame(matrix(c("Age (Years)",paste0(round(mean(NHANES_DM2_Diab$AGE),2), " (" ,"±", 
                                                            round(sd(NHANES_DM2_Diab$AGE),2),")")),ncol=2,byrow=T)),columns)

## Ethnicity 
etnicity <- table(NHANES_DM2_Diab$RACE_REC)[c(1:3)]
etnicity.prop <- round(prop.table(table(NHANES_DM2_Diab$RACE_REC)),4)[c(1:3)]*100
Race<-`names<-`(as.data.frame(matrix(c(c("Mexican American (%)","Non-Hispanic White (%)","Non-Hispanic Black (%)"),
                                       paste(etnicity,paste0('(',etnicity.prop,')'))),ncol=2,byrow=F)),columns)

#Education Categories

edu <- table(NHANES_DM2_Diab$EDU_REC_2)[c(1:3)]
edu.prop <- round(prop.table(table(NHANES_DM2_Diab$EDU_REC_2)),4)[c(1:3)]*100
Education<-`names<-`(as.data.frame(matrix(c(c("Primary School \nor \nNo-Education (%)",
                                              "Secondary, High-School \nor \nAA Degreee (%)",
                                              "College \nor \nHigher (%)"),
                                            paste(edu,paste0('(',edu.prop,')'))),ncol=2,byrow=F)),columns)


#Glycated hemoglobin

A1C<-`names<-`(as.data.frame(matrix(c("Glycated Hemoglobin (%)",paste0(round(mean(NHANES_DM2_Diab$Glycohemoglobin,na.rm=T),2), " (" ,"±", 
                                                                       round(sd(NHANES_DM2_Diab$Glycohemoglobin,na.rm=T),2),")")),ncol=2,byrow=T)),columns)


#Glucose
nortest::ad.test(NHANES_DM2_Diab$Glucose)
Glucose<-`names<-`(as.data.frame(matrix(c("Glucose (mg/dl)",paste0(round(mean(NHANES_DM2_Diab$Glucose,na.rm=T),2), " (" ,"±", 
                                                                   round(sd(NHANES_DM2_Diab$Glucose,na.rm=T),2),")")),ncol=2,byrow=T)),columns)

#Insulin
nortest::ad.test(NHANES_DM2_Diab$Insulin)

Insulin<-`names<-`(as.data.frame(matrix(c("Insulin (mUI)",paste0(round(mean(NHANES_DM2_Diab$Insulin,na.rm=T),2), " (" ,"±", 
                                                                 round(sd(NHANES_DM2_Diab$Insulin,na.rm=T),2),")")),ncol=2,byrow=T)),columns)

#Weight

Weight<-`names<-`(as.data.frame(matrix(c("Weight (Kg)",paste0(round(mean(NHANES_DM2_Diab$Weight,na.rm=T),2), " (" ,"±", 
                                                              round(sd(NHANES_DM2_Diab$Weight,na.rm=T),2),")")),ncol=2,byrow=T)),columns)

#Height
nortest::ad.test(NHANES_DM2_Diab$Height)

Height<-`names<-`(as.data.frame(matrix(c("Height (cm)",paste0(round(mean(NHANES_DM2_Diab$Height,na.rm=T),2), " (" ,"±", 
                                                              round(sd(NHANES_DM2_Diab$Height,na.rm=T),2),")")),ncol=2,byrow=T)),columns)

#BMI
nortest::ad.test(NHANES_DM2_Diab$BMI_2)

Body_Mass_Index<-`names<-`(as.data.frame(matrix(c("Body Mass Index (kg/m2)",paste0(round(mean(NHANES_DM2_Diab$BMI_2, na.rm=T),2), " (" ,"±", 
                                                                                   round(sd(NHANES_DM2_Diab$BMI_2, na.rm=T),2),")")),ncol=2,byrow=T)),columns)

#HOMA-B
nortest::ad.test(NHANES_DM2_Diab$HOMA2_B)
HOMA_B<-`names<-`(as.data.frame(matrix(c("HOMA2-B",paste0(round(mean(NHANES_DM2_Diab$HOMA2_B,na.rm=T),2), " (" ,"±", 
                                                          round(sd(NHANES_DM2_Diab$HOMA2_B,na.rm=T),2),")")),ncol=2,byrow=T)),columns)


#HOMA2-IR
nortest::ad.test(NHANES_DM2_Diab$HOMA2_IR)
HOMA_IR<-`names<-`(as.data.frame(matrix(c("HOMA2-IR",paste0(round(mean(NHANES_DM2_Diab$HOMA2_IR,na.rm=T),2), " (" ,"±", 
                                                            round(sd(NHANES_DM2_Diab$HOMA2_IR,na.rm=T),2),")")),ncol=2,byrow=T)),columns)

## Diabetes 

diab <- table(NHANES_DM2_Diab$DM2_final)[c(1)]
diab.prop <- round(prop.table(table(NHANES_DM2_Diab$DM2_final)),4)[c(1)]*100
Diabetes<-`names<-`(as.data.frame(matrix(c("Diabetes (%)",paste(diab,paste0('(',diab.prop,')'))),ncol=2,byrow=T)),columns)

## New-Onset Diabetes

new.onset <- table(NHANES_DM2_Diab$time_diabetes<5)[c(2)]
new.onset.prop <- round(prop.table(table(NHANES_DM2_Diab$time_diabetes<5)),4)[c(2)]*100
new.onset_diabetes<-`names<-`(as.data.frame(matrix(c("Recently Diagnosed Diabetes (%)",paste(new.onset,paste0('(',new.onset.prop,')'))),ncol=2,byrow=T)),columns)

## Previous Diabetes Diagnosis 

diab.diag <- table(NHANES_DM2_Diab$Dr_diabetes==1)[c(1)]
diab.diag.prop <- round(prop.table(table(NHANES_DM2_Diab$Dr_diabetes==1)),4)[c(1)]*100
Diabetes_diagnosis<-`names<-`(as.data.frame(matrix(c("Previous Diabetes Diagnosis (%)",paste(diab.diag,paste0('(',diab.diag.prop,')'))),ncol=2,byrow=T)),columns)

#Age diabetes

nortest::ad.test(NHANES_DM2_Diab$Age_diabetes)
NHANES_DM2_Diab$time_diabetes<-(NHANES_DM2_Diab$AGE-NHANES_DM2_Diab$Age_diabetes)

Diabetes_time<-`names<-`(as.data.frame(matrix(c("Time Since Diabetes Diagnosis (Years)",paste0(round(mean(NHANES_DM2_Diab$time_diabetes,na.rm=T),2), " (" ,"±", 
                                                                                               round(sd(NHANES_DM2_Diab$time_diabetes,na.rm=T),2),")")),ncol=2,byrow=T)),columns)
n  
## Diabetes Subgroup 

clus <- table(NHANES_DM2_Diab$Cluster)[c(1:4)]
clus.prop <- round(prop.table(table(NHANES_DM2_Diab$Cluster)),4)[c(1:4)]*100

Cluster<-`names<-`(as.data.frame(matrix(c(c("MARD","MOD","SIID","SIRD"),
                                          paste(clus,paste0('(',clus.prop,')'))),ncol=2,byrow=F)),columns)

#1.3 Table 1 Sample Characteristics

Table1<-rbind(Sexo,Edad,Race,Education,A1C,Glucose,Insulin,Weight,Height,Body_Mass_Index,HOMA_B,HOMA_IR, Diabetes_time,new.onset_diabetes, Cluster)
Table1<-flextable::align(flextable::flextable(Table1,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table1,path="Table_1.docx")


#----Stratification of NHAHES cycles (diabetes and clycle) (Supplementary Table 3)----

NHANES_DM2_Diab<-NHANES_DM2%>%dplyr::filter(AGE_18==1)%>%dplyr::filter(DM2_final==1)%>%dplyr::filter(!is.na(Cluster))
columns <- c('Parameter',"1988-1990 (1020)","1991-1994 (689)", "1999-2002 (1280)", "2003-2006 (1311)",
             "2007-2010 (1991)", "2011-2014 (1852)","2015-2018 (2172)")

## Sex Female

sexo <- table(NHANES_DM2_Diab$GENDER_REC,NHANES_DM2_Diab$YEAR_2,useNA = "always")[c(1,4,7,10,13,16,19)]
sexoprop <- round(prop.table(table(NHANES_DM2_Diab$GENDER_REC,NHANES_DM2_Diab$YEAR_2,useNA = "always"),2),4)[c(1,4,7,10,13,16,19)]*100

Sexo<-`names<-`(as.data.frame(matrix(c("Female (%)",paste(sexo,paste0('(',sexoprop,')'))),ncol=8,byrow=T)),columns)

#Age

Edad<-`names<-`(as.data.frame(matrix(c("Age (Years)",paste0(round(tapply(NHANES_DM2_Diab$AGE,NHANES_DM2_Diab$YEAR_2, mean),2), " (" ,"±", 
                                                           round(tapply(NHANES_DM2_Diab$AGE,NHANES_DM2_Diab$YEAR_2, sd),2),")")),ncol=8,byrow=T)),columns)

## Ethnicity 
etnicity <- table(NHANES_DM2_Diab$RACE_REC,NHANES_DM2_Diab$YEAR_2)[c(1:21)]
etnicity.prop <- round(prop.table(table(NHANES_DM2_Diab$RACE_REC,NHANES_DM2_Diab$YEAR_2),2),4)[c(1:21)]*100
Race<-`names<-`(as.data.frame(matrix(c(c("Mexican American (%)","Non-Hispanic White (%)","Non-Hispanic Black (%)"),
                                            paste(etnicity,paste0('(',etnicity.prop,')'))),ncol=8,byrow=F)),columns)

#Education Categories

edu <- table(NHANES_DM2_Diab$EDU_REC_2,NHANES_DM2_Diab$YEAR_2)[c(1:21)]
edu.prop <- round(prop.table(table(NHANES_DM2_Diab$EDU_REC_2,NHANES_DM2_Diab$YEAR_2),2),4)[c(1:21)]*100
Education<-`names<-`(as.data.frame(matrix(c(c("Primary School \nor \nNo-Education (%)",
                                              "Secondary, High-School \nor \nAA Degreee (%)",
                                              "College \nor \nHigher (%)"),
                                            paste(edu,paste0('(',edu.prop,')'))),ncol=8,byrow=F)),columns)


#Glycated hemoglobin

A1C<-`names<-`(as.data.frame(matrix(c("Glycated Hemoglobin (%)",paste0(round(tapply(NHANES_DM2_Diab$Glycohemoglobin,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                           round(tapply(NHANES_DM2_Diab$Glycohemoglobin,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)


#Glucose
nortest::ad.test(NHANES_DM2_Diab$Glucose)

Glucose<-`names<-`(as.data.frame(matrix(c("Glucose (mg/dl)",paste0(round(tapply(NHANES_DM2_Diab$Glycohemoglobin,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                               round(tapply(NHANES_DM2_Diab$Glycohemoglobin,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)

#Insulin
nortest::ad.test(NHANES_DM2_Diab$Insulin)

Insulin<-`names<-`(as.data.frame(matrix(c("Insulin (mUI)",paste0(round(tapply(NHANES_DM2_Diab$Insulin,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                                   round(tapply(NHANES_DM2_Diab$Insulin,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)

#Weight

Weight<-`names<-`(as.data.frame(matrix(c("Weight (Kg)",paste0(round(tapply(NHANES_DM2_Diab$Weight,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                                 round(tapply(NHANES_DM2_Diab$Weight,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)


#Height
nortest::ad.test(NHANES_DM2_Diab$Height)

Height<-`names<-`(as.data.frame(matrix(c("Height (cm)",paste0(round(tapply(NHANES_DM2_Diab$Height,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab$Height,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)

#BMI
nortest::ad.test(NHANES_DM2_Diab$BMI_2)

Body_Mass_Index<-`names<-`(as.data.frame(matrix(c("Body Mass Index (kg/m2)",paste0(round(tapply(NHANES_DM2_Diab$BMI_2,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab$BMI_2,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)

#HOMA-B
nortest::ad.test(NHANES_DM2_Diab$HOMA2_B)
HOMA_B<-`names<-`(as.data.frame(matrix(c("HOMA2-B",paste0(round(tapply(NHANES_DM2_Diab$HOMA2_B,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                                                   round(tapply(NHANES_DM2_Diab$HOMA2_B,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)


#HOMA2-IR
nortest::ad.test(NHANES_DM2_Diab$HOMA2_IR)
HOMA_IR<-`names<-`(as.data.frame(matrix(c("HOMA2-IR",paste0(round(tapply(NHANES_DM2_Diab$HOMA2_IR,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                          round(tapply(NHANES_DM2_Diab$HOMA2_IR,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)

## Diabetes 

diab <- table(NHANES_DM2_Diab$DM2_final,NHANES_DM2_Diab$YEAR_2)[c(1:7)]
diab.prop <- round(prop.table(table(NHANES_DM2_Diab$DM2_final,NHANES_DM2_Diab$YEAR_2),2),4)[c(1:7)]*100
Diabetes<-`names<-`(as.data.frame(matrix(c("Diabetes (%)",paste(diab,paste0('(',diab.prop,')'))),ncol=8,byrow=T)),columns)

## New-Onset Diabetes

new.onset <- table(NHANES_DM2_Diab$time_diabetes<5,NHANES_DM2_Diab$YEAR_2)[c(2,4,6,8,10,12,14)]
new.onset.prop <- round(prop.table(table(NHANES_DM2_Diab$time_diabetes<5,NHANES_DM2_Diab$YEAR_2),2),4)[c(2,4,6,8,10,12,14)]*100

new.onset_diabetes<-`names<-`(as.data.frame(matrix(c("Recently Diagnosed Diabetes (%)",paste(new.onset,paste0('(',new.onset.prop,')'))),ncol=8,byrow=T)),columns)

## Previous Diabetes Diagnosis 

diab.diag <- table(NHANES_DM2_Diab$Dr_diabetes==1,NHANES_DM2_Diab$YEAR_2)[c(1:7)]
diab.diag.prop <- round(prop.table(table(NHANES_DM2_Diab$Dr_diabetes==1,NHANES_DM2_Diab$YEAR_2),2),4)[c(1:7)]*100
Diabetes_diagnosis<-`names<-`(as.data.frame(matrix(c("Previous Diabetes Diagnosis (%)",paste(diab.diag,paste0('(',diab.diag.prop,')'))),ncol=8,byrow=T)),columns)

#Age diabetes

nortest::ad.test(NHANES_DM2_Diab$Age_diabetes)
NHANES_DM2_Diab$time_diabetes<-(NHANES_DM2_Diab$AGE-NHANES_DM2_Diab$Age_diabetes)

Diabetes_time<-`names<-`(as.data.frame(matrix(c("Time Since Diabetes Diagnosis (Years)",paste0(round(tapply(NHANES_DM2_Diab$time_diabetes,NHANES_DM2_Diab$YEAR_2, mean,na.rm=T),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab$time_diabetes,NHANES_DM2_Diab$YEAR_2, sd,na.rm=T),2),")")),ncol=8,byrow=T)),columns)

## Diabetes Subgroup 

clus <- table(NHANES_DM2_Diab$Cluster,NHANES_DM2_Diab$YEAR_2)[c(1:28)]
clus.prop <- round(prop.table(table(NHANES_DM2_Diab$Cluster,NHANES_DM2_Diab$YEAR_2),2),4)[c(1:28)]*100

Cluster<-`names<-`(as.data.frame(matrix(c(c("MARD","MOD","SIID","SIRD"),
                                            paste(clus,paste0('(',clus.prop,')'))),ncol=8,byrow=F)),columns)

#1.3 Table 1 Sample Characteristics

Table3<-rbind(Sexo,Edad,Race,Education,A1C,Glucose,Insulin,Weight,Height,Body_Mass_Index,HOMA_B,HOMA_IR, Diabetes_time,new.onset_diabetes, Cluster)
Table3<-flextable::align(flextable::flextable(Table3,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table3,path="Supplementary_Table_3.docx")


#----Stratification by phenotypes (diabetes and cluster) (Supplementary Table 4)----

NHANES_DM2_Diab_Clus<-NHANES_DM2%>%dplyr::filter(AGE_18==1)%>%dplyr::filter(DM2_final==1)%>%dplyr::filter(!is.na(Cluster))
columns <- c('Parameter',"MARD (1403)","MOD (1654)", "SIID (1290)", "SIRD (1142)")

## Sex Female

sexo <- table(NHANES_DM2_Diab_Clus$GENDER_REC,NHANES_DM2_Diab_Clus$Cluster)[c(1,3,5,7)]
sexoprop <- round(prop.table(table(NHANES_DM2_Diab_Clus$GENDER_REC,NHANES_DM2_Diab_Clus$Cluster),2),4)[c(1,3,5,7)]*100

Sexo<-`names<-`(as.data.frame(matrix(c("Female (%)",paste(sexo,paste0('(',sexoprop,')'))),ncol=5,byrow=T)),columns)

#Age

Edad<-`names<-`(as.data.frame(matrix(c("Age (Years)",paste0(round(tapply(NHANES_DM2_Diab_Clus$AGE,NHANES_DM2_Diab_Clus$Cluster, mean),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_Clus$AGE,NHANES_DM2_Diab_Clus$Cluster, sd),2),")")),ncol=5,byrow=T)),columns)

## Ethnicity 
etnicity <- table(NHANES_DM2_Diab_Clus$RACE_REC,NHANES_DM2_Diab_Clus$Cluster)[c(1:12)]
etnicity.prop <- round(prop.table(table(NHANES_DM2_Diab_Clus$RACE_REC,NHANES_DM2_Diab_Clus$YEAR_2),2),4)[c(1:12)]*100
Race<-`names<-`(as.data.frame(matrix(c(c("Mexican American (%)","Non-Hispanic White (%)","Non-Hispanic Black (%)"),
                                       paste(etnicity,paste0('(',etnicity.prop,')'))),ncol=5,byrow=F)),columns)

#Education Categories

edu <- table(NHANES_DM2_Diab_Clus$EDU_REC_2,NHANES_DM2_Diab_Clus$Cluster)[c(1:12)]
edu.prop <- round(prop.table(table(NHANES_DM2_Diab_Clus$EDU_REC_2,NHANES_DM2_Diab_Clus$Cluster),2),4)[c(1:12)]*100
Education<-`names<-`(as.data.frame(matrix(c(c("Primary School \nor \nNo-Education Secondary(%)",
                                              "High-School \nor \nAA Degreee (%)",
                                              "College \nor \nHigher (%)"),
                                            paste(edu,paste0('(',edu.prop,')'))),ncol=5,byrow=F)),columns)


#Glycated hemoglobin

A1C<-`names<-`(as.data.frame(matrix(c("Glycated Hemoglobin (%)",paste0(round(tapply(NHANES_DM2_Diab_Clus$Glycohemoglobin,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                                       round(tapply(NHANES_DM2_Diab_Clus$Glycohemoglobin,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)


#Glucose
nortest::ad.test(NHANES_DM2_Diab_Clus$Glucose)

Glucose<-`names<-`(as.data.frame(matrix(c("Glucose (mg/dl)",paste0(round(tapply(NHANES_DM2_Diab_Clus$Glucose,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                                   round(tapply(NHANES_DM2_Diab_Clus$Glucose,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)

#Insulin
nortest::ad.test(NHANES_DM2_Diab_Clus$Insulin)

Insulin<-`names<-`(as.data.frame(matrix(c("Insulin (mUI)",paste0(round(tapply(NHANES_DM2_Diab_Clus$Insulin,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                                 round(tapply(NHANES_DM2_Diab_Clus$Insulin,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)

#Weight

Weight<-`names<-`(as.data.frame(matrix(c("Weight (Kg)",paste0(round(tapply(NHANES_DM2_Diab_Clus$Weight,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_Clus$Weight,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)


#Height
nortest::ad.test(NHANES_DM2_Diab_Clus$Height)

Height<-`names<-`(as.data.frame(matrix(c("Height (cm)",paste0(round(tapply(NHANES_DM2_Diab_Clus$Height,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_Clus$Height,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)

#BMI
nortest::ad.test(NHANES_DM2_Diab_Clus$BMI_2)

Body_Mass_Index<-`names<-`(as.data.frame(matrix(c("Body Mass Index (kg/m2)",paste0(round(tapply(NHANES_DM2_Diab_Clus$BMI_2,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                                                   round(tapply(NHANES_DM2_Diab_Clus$BMI_2,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)

#HOMA-B
nortest::ad.test(NHANES_DM2_Diab_Clus$HOMA2_B)
HOMA_B<-`names<-`(as.data.frame(matrix(c("HOMA2-B",paste0(round(tapply(NHANES_DM2_Diab_Clus$HOMA2_B,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                          round(tapply(NHANES_DM2_Diab_Clus$HOMA2_B,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)


#HOMA2-IR
nortest::ad.test(NHANES_DM2_Diab_Clus$HOMA2_IR)
HOMA_IR<-`names<-`(as.data.frame(matrix(c("HOMA2-IR",paste0(round(tapply(NHANES_DM2_Diab_Clus$HOMA2_IR,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_Clus$HOMA2_IR,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)


## Previous Diabetes Diagnosis 

diab.diag <- table(NHANES_DM2_Diab_Clus$Dr_diabetes==1,NHANES_DM2_Diab_Clus$Cluster)[c(2,4,6,8)]
diab.diag.prop <- round(prop.table(table(NHANES_DM2_Diab_Clus$Dr_diabetes==1,NHANES_DM2_Diab_Clus$Cluster),2),4)[c(2,4,6,8)]*100
Diabetes_diagnosis<-`names<-`(as.data.frame(matrix(c("Previous Diabetes Diagnosis (%)",paste(diab.diag,paste0('(',diab.diag.prop,')'))),ncol=5,byrow=T)),columns)

## New-Onset Diabetes

new.onset <- table(NHANES_DM2_Diab_Clus$time_diabetes<5,NHANES_DM2_Diab_Clus$Cluster)[c(2,4,6,8)]
new.onset.prop <- round(prop.table(table(NHANES_DM2_Diab$time_diabetes<5,NHANES_DM2_Diab$Cluster),2),4)[c(2,4,6,8)]*100
new.onset_diabetes<-`names<-`(as.data.frame(matrix(c("Diabetes (%)",paste(new.onset,paste0('(',new.onset.prop,')'))),ncol=5,byrow=T)),columns)


#Age diabetes

nortest::ad.test(NHANES_DM2_Diab_Clus$Age_diabetes)
NHANES_DM2_Diab_Clus$time_diabetes<-(NHANES_DM2_Diab_Clus$AGE-NHANES_DM2_Diab_Clus$Age_diabetes)

Diabetes_time<-`names<-`(as.data.frame(matrix(c("Time Since Diabetes Diagnosis (Years)",paste0(round(tapply(NHANES_DM2_Diab_Clus$time_diabetes,NHANES_DM2_Diab_Clus$Cluster, mean,na.rm=T),2), " (" ,"±", 
                                                                                               round(tapply(NHANES_DM2_Diab_Clus$time_diabetes,NHANES_DM2_Diab_Clus$Cluster, sd,na.rm=T),2),")")),ncol=5,byrow=T)),columns)

#1.3 Table 1 Sample Characteristics

Table4<-rbind(Sexo,Edad,Race,Education,A1C,Glucose,Insulin,Weight,Height,Body_Mass_Index,HOMA_B,HOMA_IR,Diabetes_diagnosis, Diabetes_time)
Table4<-flextable::align(flextable::flextable(Table4,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table4,path="Supplementary_Table_4.docx")


#----Stratification by ethinicity (diabetes and ethnicity) (Supplementary Table 5)----

NHANES_DM2_Diab_Race<-NHANES_DM2%>%dplyr::filter(AGE_18==1)%>%dplyr::filter(DM2_final==1)%>%dplyr::filter(!is.na(RACE_REC))
columns <- c('Parameter',"Mexican American  (3008)","Non-Hispanic White (4516)", "Non-Hispanic Black (2791)")

## Sex Female

sexo <- table(NHANES_DM2_Diab_Race$GENDER_REC,NHANES_DM2_Diab_Race$RACE_REC)[c(1,3,5)]
sexoprop <- round(prop.table(table(NHANES_DM2_Diab_Race$GENDER_REC,NHANES_DM2_Diab_Race$RACE_REC),2),4)[c(1,3,5)]*100
Sexo<-`names<-`(as.data.frame(matrix(c("Female (%)",paste(sexo,paste0('(',sexoprop,')'))),ncol=4,byrow=T)),columns)

#Age
Edad<-`names<-`(as.data.frame(matrix(c("Age (Years)",paste0(round(tapply(NHANES_DM2_Diab_Race$AGE,NHANES_DM2_Diab_Race$RACE_REC, mean),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_Race$AGE,NHANES_DM2_Diab_Race$RACE_REC, sd),2),")")),ncol=4,byrow=T)),columns)

#Education Categories

edu <- table(NHANES_DM2_Diab_Race$EDU_REC_2,NHANES_DM2_Diab_Race$RACE_REC)[c(1:9)]
edu.prop <- round(prop.table(table(NHANES_DM2_Diab_Race$EDU_REC_2,NHANES_DM2_Diab_Race$RACE_REC),2),4)[c(1:9)]*100
Education<-`names<-`(as.data.frame(matrix(c(c("Primary School \nor \nNo-Education Secondary(%)",
                                              "High-School \nor \nAA Degreee (%)",
                                              "College \nor \nHigher (%)"),
                                            paste(edu,paste0('(',edu.prop,')'))),ncol=4,byrow=F)),columns)

#Glycated hemoglobin
A1C<-`names<-`(as.data.frame(matrix(c("Glycated Hemoglobin (%)",paste0(round(tapply(NHANES_DM2_Diab_Race$Glycohemoglobin,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                       round(tapply(NHANES_DM2_Diab_Race$Glycohemoglobin,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


#Glucose
nortest::ad.test(NHANES_DM2_Diab_Race$Glucose)

Glucose<-`names<-`(as.data.frame(matrix(c("Glucose (mg/dl)",paste0(round(tapply(NHANES_DM2_Diab_Race$Glucose,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                   round(tapply(NHANES_DM2_Diab_Race$Glucose,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#Insulin
nortest::ad.test(NHANES_DM2_Diab_Race$Insulin)

Insulin<-`names<-`(as.data.frame(matrix(c("Insulin (mUI)",paste0(round(tapply(NHANES_DM2_Diab_Race$Insulin,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                 round(tapply(NHANES_DM2_Diab_Race$Insulin,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#Weight

Weight<-`names<-`(as.data.frame(matrix(c("Weight (Kg)",paste0(round(tapply(NHANES_DM2_Diab_Race$Weight,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_Race$Weight,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


#Height
nortest::ad.test(NHANES_DM2_Diab_Race$Height)

Height<-`names<-`(as.data.frame(matrix(c("Height (cm)",paste0(round(tapply(NHANES_DM2_Diab_Race$Height,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_Race$Height,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#BMI
nortest::ad.test(NHANES_DM2_Diab_Race$BMI_2)

Body_Mass_Index<-`names<-`(as.data.frame(matrix(c("Body Mass Index (kg/m2)",paste0(round(tapply(NHANES_DM2_Diab_Race$BMI_2,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                                   round(tapply(NHANES_DM2_Diab_Race$BMI_2,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#HOMA-B
nortest::ad.test(NHANES_DM2_Diab_Race$HOMA2_B)
HOMA_B<-`names<-`(as.data.frame(matrix(c("HOMA2-B",paste0(round(tapply(NHANES_DM2_Diab_Race$HOMA2_B,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                          round(tapply(NHANES_DM2_Diab_Race$HOMA2_B,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


#HOMA2-IR
nortest::ad.test(NHANES_DM2_Diab_Race$HOMA2_IR)
HOMA_IR<-`names<-`(as.data.frame(matrix(c("HOMA2-IR",paste0(round(tapply(NHANES_DM2_Diab_Race$HOMA2_IR,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_Race$HOMA2_IR,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


## Previous Diabetes Diagnosis 

diab.diag <- table(NHANES_DM2_Diab_Race$Dr_diabetes==1,NHANES_DM2_Diab_Race$RACE_REC)[c(2,4,6)]
diab.diag.prop <- round(prop.table(table(NHANES_DM2_Diab_Race$Dr_diabetes==1,NHANES_DM2_Diab_Race$RACE_REC),2),4)[c(2,4,6)]*100
Diabetes_diagnosis<-`names<-`(as.data.frame(matrix(c("Previous Diabetes Diagnosis (%)",paste(diab.diag,paste0('(',diab.diag.prop,')'))),ncol=4,byrow=T)),columns)

## New-Onset Diabetes
chisq.test(table(NHANES_DM2_Diab_Race$time_diabetes<5,NHANES_DM2_Diab_Race$RACE_REC))
new.onset <- table(NHANES_DM2_Diab_Race$time_diabetes<5,NHANES_DM2_Diab_Race$RACE_REC)[c(2,4,6)]
new.onset.prop <- round(prop.table(table(NHANES_DM2_Diab_Race$time_diabetes<5,NHANES_DM2_Diab_Race$RACE_REC),2),4)[c(2,4,6)]*100
new.onset_diabetes<-`names<-`(as.data.frame(matrix(c("Diabetes (%)",paste(new.onset,paste0('(',new.onset.prop,')'))),ncol=4,byrow=T)),columns)


#Age diabetes

nortest::ad.test(NHANES_DM2_Diab_Race$Age_diabetes)
NHANES_DM2_Diab_Race$time_diabetes<-(NHANES_DM2_Diab_Race$AGE-NHANES_DM2_Diab_Race$Age_diabetes)

Diabetes_time<-`names<-`(as.data.frame(matrix(c("Time Since Diabetes Diagnosis (Years)",paste0(round(tapply(NHANES_DM2_Diab_Race$time_diabetes,NHANES_DM2_Diab_Race$RACE_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                                               round(tapply(NHANES_DM2_Diab_Race$time_diabetes,NHANES_DM2_Diab_Race$RACE_REC, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#Diabetes Subgroup
clus <- table(NHANES_DM2_Diab_Race$Cluster,NHANES_DM2_Diab_Race$RACE_REC)[c(1:12)]
clus.prop <- round(prop.table(table(NHANES_DM2_Diab_Race$Cluster,NHANES_DM2_Diab_Race$RACE_REC),2),4)[c(1:12)]*100
Cluster<-`names<-`(as.data.frame(matrix(c(c("MARD","MOD","SIID","SIRD"),
                                            paste(clus,paste0('(',clus.prop,')'))),ncol=4,byrow=F)),columns)

#1.3 Table 1 Sample Characteristics

Table5<-rbind(Sexo,Edad,Education,A1C,Glucose,Insulin,Weight,Height,Body_Mass_Index,HOMA_B,HOMA_IR,Diabetes_diagnosis, Diabetes_time,Cluster)
Table5<-flextable::align(flextable::flextable(Table5,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table5,path="Supplementary_Table_5.docx")

#----Stratification by sex (diabetes and sex) (Supplementary Table 6)-----

NHANES_DM2_Diab_Sex<-NHANES_DM2%>%dplyr::filter(AGE_18==1)%>%dplyr::filter(DM2_final==1)%>%dplyr::filter(!is.na(GENDER_REC))

columns <- c('Parameter',"Women (5055)","Men (5260)")

#Age
Edad<-`names<-`(as.data.frame(matrix(c("Age (Years)",paste0(round(tapply(NHANES_DM2_Diab_Sex$AGE,NHANES_DM2_Diab_Sex$GENDER_REC, mean),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_Sex$AGE,NHANES_DM2_Diab_Sex$GENDER_REC, sd),2),")")),ncol=3,byrow=T)),columns)


## Ethnicity 

chisq.test(table(NHANES_DM2_Diab_Sex$RACE_REC=="Mexican American",NHANES_DM2_Diab_Sex$GENDER_REC))
chisq.test(table(NHANES_DM2_Diab_Sex$RACE_REC=="Non-Hispanic White",NHANES_DM2_Diab_Sex$GENDER_REC))
chisq.test(table(NHANES_DM2_Diab_Sex$RACE_REC=="Non-Hispanic Black",NHANES_DM2_Diab_Sex$GENDER_REC))

etnicity <- table(NHANES_DM2_Diab_Sex$RACE_REC,NHANES_DM2_Diab_Sex$GENDER_REC)[c(1:6)]
etnicity.prop <- round(prop.table(table(NHANES_DM2_Diab_Sex$RACE_REC,NHANES_DM2_Diab_Sex$GENDER_REC),2),4)[c(1:6)]*100
Race<-`names<-`(as.data.frame(matrix(c(c("Mexican American (%)","Non-Hispanic White (%)","Non-Hispanic Black (%)"),
                                       paste(etnicity,paste0('(',etnicity.prop,')'))),ncol=3,byrow=F)),columns)
#Education Categories

edu <- table(NHANES_DM2_Diab_Sex$EDU_REC_2,NHANES_DM2_Diab_Sex$GENDER_REC)[c(1:6)]
edu.prop <- round(prop.table(table(NHANES_DM2_Diab_Sex$EDU_REC_2,NHANES_DM2_Diab_Sex$GENDER_REC),2),4)[c(1:6)]*100
Education<-`names<-`(as.data.frame(matrix(c(c("Primary School \nor \nNo-Education Secondary(%)",
                                              "High-School \nor \nAA Degreee (%)",
                                              "College \nor \nHigher (%)"),
                                            paste(edu,paste0('(',edu.prop,')'))),ncol=3,byrow=F)),columns)

#Glycated hemoglobin
A1C<-`names<-`(as.data.frame(matrix(c("Glycated Hemoglobin (%)",paste0(round(tapply(NHANES_DM2_Diab_Sex$Glycohemoglobin,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                       round(tapply(NHANES_DM2_Diab_Sex$Glycohemoglobin,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)


#Glucose
nortest::ad.test(NHANES_DM2_Diab_Sex$Glucose)

Glucose<-`names<-`(as.data.frame(matrix(c("Glucose (mg/dl)",paste0(round(tapply(NHANES_DM2_Diab_Sex$Glucose,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                   round(tapply(NHANES_DM2_Diab_Sex$Glucose,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)

#Insulin

nortest::ad.test(NHANES_DM2_Diab_Sex$Insulin)

Insulin<-`names<-`(as.data.frame(matrix(c("Insulin (mUI)",paste0(round(tapply(NHANES_DM2_Diab_Sex$Insulin,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                 round(tapply(NHANES_DM2_Diab_Sex$Insulin,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)

#Weight

Weight<-`names<-`(as.data.frame(matrix(c("Weight (Kg)",paste0(round(tapply(NHANES_DM2_Diab_Sex$Weight,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_Sex$Weight,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)


#Height
nortest::ad.test(NHANES_DM2_Diab_Sex$Height)

Height<-`names<-`(as.data.frame(matrix(c("Height (cm)",paste0(round(tapply(NHANES_DM2_Diab_Sex$Height,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_Sex$Height,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)

#BMI
nortest::ad.test(NHANES_DM2_Diab_Sex$BMI_2)

Body_Mass_Index<-`names<-`(as.data.frame(matrix(c("Body Mass Index (kg/m2)",paste0(round(tapply(NHANES_DM2_Diab_Sex$BMI_2,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                                   round(tapply(NHANES_DM2_Diab_Sex$BMI_2,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)

#HOMA-B
nortest::ad.test(NHANES_DM2_Diab_Sex$HOMA2_B)
HOMA_B<-`names<-`(as.data.frame(matrix(c("HOMA2-B",paste0(round(tapply(NHANES_DM2_Diab_Sex$HOMA2_B,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                          round(tapply(NHANES_DM2_Diab_Sex$HOMA2_B,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)


#HOMA2-IR
nortest::ad.test(NHANES_DM2_Diab_Sex$HOMA2_IR)
HOMA_IR<-`names<-`(as.data.frame(matrix(c("HOMA2-IR",paste0(round(tapply(NHANES_DM2_Diab_Sex$HOMA2_IR,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_Sex$HOMA2_IR,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)


## Previous Diabetes Diagnosis 

diab.diag <- table(NHANES_DM2_Diab_Sex$Dr_diabetes==1,NHANES_DM2_Diab_Sex$GENDER_REC)[c(2,4)]
diab.diag.prop <- round(prop.table(table(NHANES_DM2_Diab_Sex$Dr_diabetes==1,NHANES_DM2_Diab_Sex$GENDER_REC),2),4)[c(2,4)]*100
Diabetes_diagnosis<-`names<-`(as.data.frame(matrix(c("Previous Diabetes Diagnosis (%)",paste(diab.diag,paste0('(',diab.diag.prop,')'))),ncol=3,byrow=T)),columns)

#Age diabetes

nortest::ad.test(NHANES_DM2_Diab_Sex$Age_diabetes)
NHANES_DM2_Diab_Sex$time_diabetes<-(NHANES_DM2_Diab_Sex$AGE-NHANES_DM2_Diab_Sex$Age_diabetes)

Diabetes_time<-`names<-`(as.data.frame(matrix(c("Time Since Diabetes Diagnosis (Years)",paste0(round(tapply(NHANES_DM2_Diab_Sex$time_diabetes,NHANES_DM2_Diab_Sex$GENDER_REC, mean,na.rm=T),2), " (" ,"±", 
                                                                                               round(tapply(NHANES_DM2_Diab_Sex$time_diabetes,NHANES_DM2_Diab_Sex$GENDER_REC, sd,na.rm=T),2),")")),ncol=3,byrow=T)),columns)

#New Onset Diabetes 

chisq.test(table(NHANES_DM2_Diab_Sex$time_diabetes<5,NHANES_DM2_Diab_Sex$GENDER_REC))
new.onset <- table(NHANES_DM2_Diab_Sex$time_diabetes<5,NHANES_DM2_Diab_Sex$GENDER_REC)[c(2,4)]
new.onset.prop <- round(prop.table(table(NHANES_DM2_Diab_Sex$time_diabetes<5,NHANES_DM2_Diab_Sex$GENDER_REC),2),4)[c(2,4)]*100
new.onset_diabetes<-`names<-`(as.data.frame(matrix(c("Diabetes (%)",paste(new.onset,paste0('(',new.onset.prop,')'))),ncol=3,byrow=T)),columns)


#Diabetes Subgroup
clus <- table(NHANES_DM2_Diab_Sex$Cluster,NHANES_DM2_Diab_Sex$GENDER_REC)[c(1:8)]
clus.prop <- round(prop.table(table(NHANES_DM2_Diab_Sex$Cluster,NHANES_DM2_Diab_Sex$GENDER_REC),2),4)[c(1:8)]*100
Cluster<-`names<-`(as.data.frame(matrix(c(c("MARD (%)","MOD (%)","SIID (%)","SIRD (%)"),
                                          paste(clus,paste0('(',clus.prop,')'))),ncol=3,byrow=F)),columns)

#1.3 Table 1 Sample Characteristics

Table6<-rbind(Edad,Education,A1C,Glucose,Insulin,Weight,Height,Body_Mass_Index,HOMA_B,HOMA_IR,Diabetes_diagnosis, Diabetes_time,Cluster)
Table6<-flextable::align(flextable::flextable(Table6,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table6,path="Supplementary_Table_6.docx")




#----Stratification by sex (diabetes and education level) (Supplementary Table 7)-----


NHANES_DM2_Diab_EDU<-NHANES_DM2%>%dplyr::filter(AGE_18==1)%>%dplyr::filter(DM2_final==1)%>%dplyr::filter(!is.na(EDU_REC_2))
columns <- c('Parameter',"Primary School or No-Education (2508)","Secondary, High-School or AA Degree (4221)","College or Higher (3586)")

## Sex Female

sexo <- table(NHANES_DM2_Diab_EDU$GENDER_REC,NHANES_DM2_Diab_EDU$EDU_REC_2)[c(1,3,5)]
sexoprop <- round(prop.table(table(NHANES_DM2_Diab_EDU$GENDER_REC,NHANES_DM2_Diab_EDU$EDU_REC_2),2),4)[c(1,3,5)]*100
Sexo<-`names<-`(as.data.frame(matrix(c("Female (%)",paste(sexo,paste0('(',sexoprop,')'))),ncol=4,byrow=T)),columns)

#Education Level
chisq.test(table(NHANES_DM2_Diab_EDU$RACE_REC=="Mexican American",NHANES_DM2_Diab_Sex$EDU_REC_2))
chisq.test(table(NHANES_DM2_Diab_EDU$RACE_REC=="Non-Hispanic White",NHANES_DM2_Diab_Sex$EDU_REC_2))
chisq.test(table(NHANES_DM2_Diab_EDU$RACE_REC=="Non-Hispanic Black",NHANES_DM2_Diab_Sex$EDU_REC_2))

etnicity <- table(NHANES_DM2_Diab_EDU$RACE_REC,NHANES_DM2_Diab_EDU$EDU_REC_2)[c(1:9)]
etnicity.prop <- round(prop.table(table(NHANES_DM2_Diab_EDU$RACE_REC,NHANES_DM2_Diab_EDU$EDU_REC_2),2),4)[c(1:9)]*100
Race<-`names<-`(as.data.frame(matrix(c(c("Mexican American (%)","Non-Hispanic White (%)","Non-Hispanic Black (%)"),
                                       paste(etnicity,paste0('(',etnicity.prop,')'))),ncol=4,byrow=F)),columns)


#Age
Edad<-`names<-`(as.data.frame(matrix(c("Age (Years)",paste0(round(tapply(NHANES_DM2_Diab_EDU$AGE,NHANES_DM2_Diab_EDU$EDU_REC_2, mean),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_EDU$AGE,NHANES_DM2_Diab_EDU$EDU_REC_2, sd),2),")")),ncol=4,byrow=T)),columns)


#Glycated hemoglobin
A1C<-`names<-`(as.data.frame(matrix(c("Glycated Hemoglobin (%)",paste0(round(tapply(NHANES_DM2_Diab_EDU$Glycohemoglobin,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                                       round(tapply(NHANES_DM2_Diab_EDU$Glycohemoglobin,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


#Glucose
nortest::ad.test(NHANES_DM2_Diab_EDU$Glucose)

Glucose<-`names<-`(as.data.frame(matrix(c("Glucose (mg/dl)",paste0(round(tapply(NHANES_DM2_Diab_EDU$Glucose,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                                   round(tapply(NHANES_DM2_Diab_EDU$Glucose,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#Insulin

nortest::ad.test(NHANES_DM2_Diab_EDU$Insulin)

Insulin<-`names<-`(as.data.frame(matrix(c("Insulin (mUI)",paste0(round(tapply(NHANES_DM2_Diab_EDU$Insulin,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                                 round(tapply(NHANES_DM2_Diab_EDU$Insulin,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#Weight

Weight<-`names<-`(as.data.frame(matrix(c("Weight (Kg)",paste0(round(tapply(NHANES_DM2_Diab_EDU$Weight,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_EDU$Weight,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


#Height
nortest::ad.test(NHANES_DM2_Diab_EDU$Height)

Height<-`names<-`(as.data.frame(matrix(c("Height (cm)",paste0(round(tapply(NHANES_DM2_Diab_EDU$Height,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                              round(tapply(NHANES_DM2_Diab_EDU$Height,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#BMI
nortest::ad.test(NHANES_DM2_Diab_EDU$BMI_2)

Body_Mass_Index<-`names<-`(as.data.frame(matrix(c("Body Mass Index (kg/m2)",paste0(round(tapply(NHANES_DM2_Diab_EDU$BMI_2,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                                                   round(tapply(NHANES_DM2_Diab_EDU$BMI_2,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)

#HOMA-B
nortest::ad.test(NHANES_DM2_Diab_EDU$HOMA2_B)
HOMA_B<-`names<-`(as.data.frame(matrix(c("HOMA2-B",paste0(round(tapply(NHANES_DM2_Diab_EDU$HOMA2_B,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                          round(tapply(NHANES_DM2_Diab_EDU$HOMA2_B,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


#HOMA2-IR
nortest::ad.test(NHANES_DM2_Diab_EDU$HOMA2_IR)
HOMA_IR<-`names<-`(as.data.frame(matrix(c("HOMA2-IR",paste0(round(tapply(NHANES_DM2_Diab_EDU$HOMA2_IR,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                            round(tapply(NHANES_DM2_Diab_EDU$HOMA2_IR,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)


## Previous Diabetes Diagnosis 

diab.diag <- table(NHANES_DM2_Diab_EDU$Dr_diabetes==1,NHANES_DM2_Diab_EDU$EDU_REC_2)[c(2,4,6)]
diab.diag.prop <- round(prop.table(table(NHANES_DM2_Diab_EDU$Dr_diabetes==1,NHANES_DM2_Diab_EDU$EDU_REC_2),2),4)[c(2,4,6)]*100
Diabetes_diagnosis<-`names<-`(as.data.frame(matrix(c("Previous Diabetes Diagnosis (%)",paste(diab.diag,paste0('(',diab.diag.prop,')'))),ncol=4,byrow=T)),columns)

#Age diabetes

nortest::ad.test(NHANES_DM2_Diab_EDU$Age_diabetes)
NHANES_DM2_Diab_EDU$time_diabetes<-(NHANES_DM2_Diab_EDU$AGE-NHANES_DM2_Diab_EDU$Age_diabetes)
Diabetes_time<-`names<-`(as.data.frame(matrix(c("Time Since Diabetes Diagnosis (Years)",paste0(round(tapply(NHANES_DM2_Diab_EDU$time_diabetes,NHANES_DM2_Diab_EDU$EDU_REC_2, mean,na.rm=T),2), " (" ,"±", 
                                                                                               round(tapply(NHANES_DM2_Diab_EDU$time_diabetes,NHANES_DM2_Diab_EDU$EDU_REC_2, sd,na.rm=T),2),")")),ncol=4,byrow=T)),columns)
#New Onset Diabetes 

chisq.test(table(NHANES_DM2_Diab_EDU$time_diabetes<5,NHANES_DM2_Diab_EDU$EDU_REC_2))
new.onset <- table(NHANES_DM2_Diab_EDU$time_diabetes<5,NHANES_DM2_Diab_EDU$EDU_REC_2)[c(2,4,6)]
new.onset.prop <- round(prop.table(table(NHANES_DM2_Diab_EDU$time_diabetes<5,NHANES_DM2_Diab_EDU$EDU_REC_2),2),4)[c(2,4,6)]*100
new.onset_diabetes<-`names<-`(as.data.frame(matrix(c("Diabetes (%)",paste(new.onset,paste0('(',new.onset.prop,')'))),ncol=4,byrow=T)),columns)



#Diabetes Subgroup
clus <- table(NHANES_DM2_Diab_EDU$Cluster,NHANES_DM2_Diab_EDU$EDU_REC_2)[c(1:12)]
clus.prop <- round(prop.table(table(NHANES_DM2_Diab_EDU$Cluster,NHANES_DM2_Diab_EDU$EDU_REC_2),2),4)[c(1:12)]*100
Cluster<-`names<-`(as.data.frame(matrix(c(c("MARD (%)","MOD (%)","SIID (%)","SIRD (%)"),
                                          paste(clus,paste0('(',clus.prop,')'))),ncol=4,byrow=F)),columns)

#1.3 Table 1 Sample Characteristics

Table7<-rbind(Sexo,Edad,A1C,Glucose,Insulin,Weight,Height,Body_Mass_Index,HOMA_B,HOMA_IR,Diabetes_diagnosis,Diabetes_time,Cluster)
Table7<-flextable::align(flextable::flextable(Table7,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(Table7,path="Supplementary_Table_7.docx")


#----Overall diabetes trends----

DM2_prev<-svyby(~DM2_final, ~YEAR_2,NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD<-svyby(~ClusterMARD, ~YEAR_2,NHANES_subset1, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD<-svyby(~ClusterMOD, ~YEAR_2,NHANES_subset1, svymean) %>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD<-svyby(~ClusterSIRD, ~YEAR_2,NHANES_subset1, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID<-svyby(~ClusterSIID, ~YEAR_2,NHANES_subset1, svymean) %>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev)[names(DM2_prev) == "DM2_final"] <- "prop"
names(DM2_MARD)[names(DM2_MARD) == "ClusterMARD"] <- "prop"
names(DM2_MOD)[names(DM2_MOD) == "ClusterMOD"] <- "prop"
names(DM2_SIRD)[names(DM2_SIRD) == "ClusterSIRD"] <- "prop"
names(DM2_SIID)[names(DM2_SIID) == "ClusterSIID"] <- "prop"

DM2_prev$group<-"Type 2 Diabetes"
DM2_MARD$group<-"MARD"
DM2_MOD$group<-"MOD"
DM2_SIRD$group<-"SIRD"
DM2_SIID$group<-"SIDD"

prev_year<-rbind(DM2_prev,DM2_MARD,DM2_MOD,DM2_SIRD,DM2_SIID)
prev_year$subgroup<-ifelse(prev_year$group=="Type 2 Diabetes", 1, 0)
prev_year$subgroup<-factor(prev_year$subgroup, labels = c("Subgroups", "Overall"))


#----Sex stratification-----

DM2_prev_sex<-svyby(~DM2_final,~YEAR_2+GENDER_REC, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_sex<-svyby(~ClusterMARD, ~YEAR_2+GENDER_REC, NHANES_subset1, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_sex<-svyby(~ClusterMOD, ~YEAR_2+GENDER_REC , NHANES_subset1, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_sex<-svyby(~ClusterSIRD,~YEAR_2+GENDER_REC, NHANES_subset1, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_sex<-svyby(~ClusterSIID, ~YEAR_2+GENDER_REC , NHANES_subset1, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_sex)[names(DM2_prev_sex) == "DM2_final"] <- "prop"
names(DM2_MARD_sex)[names(DM2_MARD_sex) == "ClusterMARD"] <- "prop"
names(DM2_MOD_sex)[names(DM2_MOD_sex) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_sex)[names(DM2_SIRD_sex) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_sex)[names(DM2_SIID_sex) == "ClusterSIID"] <- "prop"

DM2_prev_sex$group<-"Type 2 Diabetes"
DM2_MARD_sex$group<-"MARD"
DM2_MOD_sex$group<-"MOD"
DM2_SIRD_sex$group<-"SIRD"
DM2_SIID_sex$group<-"SIDD"

prev_year_sex<-rbind(DM2_prev_sex,DM2_MARD_sex,DM2_MOD_sex,DM2_SIRD_sex,DM2_SIID_sex)
prev_year_sex$subgroup<-ifelse(prev_year_sex$group=="Type 2 Diabetes", 0, 1)
prev_year_sex$subgroup<-factor(prev_year_sex$subgroup, labels = c("Overall","Subgroups"))

prev_year_sex$low_limit<-prev_year_sex$prop-prev_year_sex$se
prev_year_sex$upp_limit<-prev_year_sex$prop+prev_year_sex$se
prev_year_sex$IC<-paste0(str_pad(prev_year_sex$prop, 1,pad = "0"),
                         c(" "),
                         c("("),
                         str_pad(prev_year_sex$low_limit, 1,pad = "0"),
                         c("-"),
                         str_pad(prev_year_sex$upp_limit,1, pad="0"),
                         c(")"))

#----Ethnicity stratification-----

DM2_prev_race<-svyby(~DM2_final,~YEAR_2+RACE_REC, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_race<-svyby(~ClusterMARD, ~YEAR_2+RACE_REC, NHANES_subset1, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_race<-svyby(~ClusterMOD, ~YEAR_2+RACE_REC , NHANES_subset1, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_race<-svyby(~ClusterSIRD,~YEAR_2+RACE_REC, NHANES_subset1, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_race<-svyby(~ClusterSIID, ~YEAR_2+RACE_REC , NHANES_subset1, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_race)[names(DM2_prev_race) == "DM2_final"] <- "prop"
names(DM2_MARD_race)[names(DM2_MARD_race) == "ClusterMARD"] <- "prop"
names(DM2_MOD_race)[names(DM2_MOD_race) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_race)[names(DM2_SIRD_race) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_race)[names(DM2_SIID_race) == "ClusterSIID"] <- "prop"

DM2_prev_race$group<-"Type 2 Diabetes"
DM2_MARD_race$group<-"MARD"
DM2_MOD_race$group<-"MOD"
DM2_SIRD_race$group<-"SIRD"
DM2_SIID_race$group<-"SIID"

prev_year_race<-rbind(DM2_prev_race,DM2_MARD_race,DM2_MOD_race,DM2_SIRD_race,DM2_SIID_race)
prev_year_race$subgroup<-ifelse(prev_year_race$group=="Type 2 Diabetes", 0, 1)
prev_year_race$subgroup<-factor(prev_year_race$subgroup, labels = c("Overall","Subgroups"))

prev_year_race$low_limit<-prev_year_race$prop-prev_year_race$se
prev_year_race$upp_limit<-prev_year_race$prop+prev_year_race$se
prev_year_race$IC<-paste0(str_pad(prev_year_race$prop, 1,pad = "0"),
                          c(" "),
                          c("("),
                          str_pad(prev_year_race$low_limit, 1,pad = "0"),
                          c("-"),
                          str_pad(prev_year_race$upp_limit,1, pad="0"),
                          c(")"))


#----Age stratification-----

DM2_prev_age<-svyby(~DM2_final,~YEAR_2+age_cat, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_age<-svyby(~ClusterMARD, ~YEAR_2+age_cat, NHANES_subset1, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_age<-svyby(~ClusterMOD, ~YEAR_2+age_cat , NHANES_subset1, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_age<-svyby(~ClusterSIRD,~YEAR_2+age_cat, NHANES_subset1, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_age<-svyby(~ClusterSIID, ~YEAR_2+age_cat , NHANES_subset1, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_age)[names(DM2_prev_age) == "DM2_final"] <- "prop"
names(DM2_MARD_age)[names(DM2_MARD_age) == "ClusterMARD"] <- "prop"
names(DM2_MOD_age)[names(DM2_MOD_age) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_age)[names(DM2_SIRD_age) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_age)[names(DM2_SIID_age) == "ClusterSIID"] <- "prop"

DM2_prev_age$group<-"Type 2 Diabetes"
DM2_MARD_age$group<-"MARD"
DM2_MOD_age$group<-"MOD"
DM2_SIRD_age$group<-"SIRD"
DM2_SIID_age$group<-"SIID"

prev_year_age<-rbind(DM2_prev_age,DM2_MARD_age,DM2_MOD_age,DM2_SIRD_age,DM2_SIID_age)
prev_year_age$subgroup<-ifelse(prev_year_age$group=="Type 2 Diabetes", 0, 1)
prev_year_age$subgroup<-factor(prev_year_age$subgroup, labels = c("Overall","Subgroups"))

prev_year_age$low_limit<-prev_year_age$prop-prev_year_age$se
prev_year_age$upp_limit<-prev_year_age$prop+prev_year_age$se
prev_year_age$IC<-paste0(str_pad(prev_year_age$prop, 1,pad = "0"),
                          c(" "),
                          c("("),
                          str_pad(prev_year_age$low_limit, 1,pad = "0"),
                          c("-"),
                          str_pad(prev_year_age$upp_limit,1, pad="0"),
                          c(")"))


#----Education stratification-----

DM2_prev_edu<-svyby(~DM2_final,~YEAR_2+EDU_REC_2, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_edu<-svyby(~ClusterMARD, ~YEAR_2+EDU_REC_2, NHANES_subset1, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_edu<-svyby(~ClusterMOD, ~YEAR_2+EDU_REC_2 , NHANES_subset1, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_edu<-svyby(~ClusterSIRD,~YEAR_2+EDU_REC_2, NHANES_subset1, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_edu<-svyby(~ClusterSIID, ~YEAR_2+EDU_REC_2 , NHANES_subset1, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_edu)[names(DM2_prev_edu) == "DM2_final"] <- "prop"
names(DM2_MARD_edu)[names(DM2_MARD_edu) == "ClusterMARD"] <- "prop"
names(DM2_MOD_edu)[names(DM2_MOD_edu) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_edu)[names(DM2_SIRD_edu) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_edu)[names(DM2_SIID_edu) == "ClusterSIID"] <- "prop"

DM2_prev_edu$group<-"Type 2 Diabetes"
DM2_MARD_edu$group<-"MARD"
DM2_MOD_edu$group<-"MOD"
DM2_SIRD_edu$group<-"SIRD"
DM2_SIID_edu$group<-"SIID"

prev_year_edu<-rbind(DM2_prev_edu,DM2_MARD_edu,DM2_MOD_edu,DM2_SIRD_edu,DM2_SIID_edu)
prev_year_edu$subgroup<-ifelse(prev_year_edu$group=="Type 2 Diabetes", 0, 1)
prev_year_edu$subgroup<-factor(prev_year_edu$subgroup, labels = c("Overall","Subgroups"))

prev_year_edu$low_limit<-prev_year_edu$prop-prev_year_edu$se
prev_year_edu$upp_limit<-prev_year_edu$prop+prev_year_edu$se
prev_year_edu$IC<-paste0(str_pad(prev_year_edu$prop, 1,pad = "0"),
                         c(" "),
                         c("("),
                         str_pad(prev_year_edu$low_limit, 1,pad = "0"),
                         c("-"),
                         str_pad(prev_year_edu$upp_limit,1, pad="0"),
                         c(")"))


#----Figure 1-----

Figure1A<-ggplot2::ggplot(prev_year[prev_year$subgroup=="Overall",], aes(x=YEAR_2, y=prop,group=1, col="Type 2 Diabetes")) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  labs(colour="")+
  scale_color_jama()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits = c(7,15),breaks = seq(7, 15, by = 1))

Figure1B<-ggplot(prev_year[prev_year$subgroup=="Subgroups",], aes(x=YEAR_2, y=prop,group=group, col=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  labs(colour="Groups")+
  scale_color_jama()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_y_continuous(limits = c(0,8),breaks = seq(0, 8, by = 1))
  geom_text(aes(x = 2,y = 8,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year[prev_year$subgroup=="Overall",]))

Figure1_FINAL<-ggpubr::ggarrange(Figure1A,Figure1B,nrow = 1, ncol = 2,common.legend = F)

ggsave(Figure1_FINAL,
       filename = "Figure1.jpg", 
       width = 25, 
       height = 10,
       units=c("cm"),
       dpi = 450,
       limitsize = FALSE)


prev_year_sex$group_2<-factor(prev_year_sex$group)
prev_year_sex$group_2<-relevel(prev_year_sex$group_2,ref = "SIID")
table(prev_year_sex$subgroup)

m1<-lme4::lmer(as.numeric(prop)~factor(GENDER_REC)*factor(group_2)+(1|YEAR_2/group_2), data=prev_year_sex[prev_year_sex$subgroup=="Subgroups",])
lmerTest::ranova(m1)
confint(m1)

Figure5<-ggpubr::ggarrange(Figure5A,Figure5_Down,nrow = 2, ncol = 1,labels = c("C"),common.legend = T)

#----Figure 2-----
#Gender Stratification
Figure2A<-ggplot(prev_year_sex[prev_year_sex$subgroup=="Overall",], aes(x=YEAR_2, y=prop, group=GENDER_REC, colour=GENDER_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Gender")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 25,label = 'p value = 0.004'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ GENDER_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_sex[prev_year_sex$subgroup=="Overall",]))

Figure2B<-ggplot(prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="MARD",], aes(x=YEAR_2, y=prop, group=GENDER_REC, colour=GENDER_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Gender")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.270'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ GENDER_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="MARD",]))

Figure2C<-ggplot(prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="MOD",], aes(x=YEAR_2, y=prop, group=GENDER_REC, colour=GENDER_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Gender")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0,10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.790'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ GENDER_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="MOD",]))

Figure2D<-ggplot(prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="SIRD",], aes(x=YEAR_2, y=prop, group=GENDER_REC, colour=GENDER_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Gender")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value = 0.220'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ GENDER_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="SIRD",]))

Figure2E<-ggplot(prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="SIDD",], aes(x=YEAR_2, y=prop, group=GENDER_REC, colour=GENDER_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Gender")+ggtitle("SIDD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value = 0.01'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ GENDER_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_sex[prev_year_sex$subgroup=="Subgroups" &prev_year_sex$group=="SIID",]))

Figure2_Down<-ggpubr::ggarrange(Figure2B,Figure2C,Figure2D,Figure2E,nrow = 2, ncol = 2,common.legend = T)
Figure2<-ggpubr::ggarrange(Figure2A,Figure2_Down,nrow = 2, ncol = 1,labels = c("A"),common.legend = T)


#Etnicity Stratificaton

Figure3A<-ggplot(prev_year_race[prev_year_race$subgroup=="Overall",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 25,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race[prev_year_race$subgroup=="Overall",]))

Figure3B<-ggplot(prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="MARD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.014'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="MARD",]))

Figure3C<-ggplot(prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="MOD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="MOD",]))

Figure3D<-ggplot(prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="SIRD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="SIRD",]))

Figure3E<-ggplot(prev_year_race[prev_year_race$subgroup=="Subgroups" & prev_year_race$group=="SIID",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIDD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ factor(RACE_REC) * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race[prev_year_race$subgroup=="Subgroups" &prev_year_race$group=="SIID",]))

Figure3_Down<-ggpubr::ggarrange(Figure3B,Figure3C,Figure3D,Figure3E,nrow = 2, ncol = 2,common.legend = T)
Figure3<-ggpubr::ggarrange(Figure3A,Figure3_Down,nrow = 2, ncol = 1,labels = c("B"),common.legend = T)

#Education Stratification

Figure5A<-ggplot(prev_year_edu[prev_year_edu$subgroup=="Overall",], 
       aes(x=YEAR_2, y=prop, group=EDU_REC_2, colour=EDU_REC_2)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Education \nLevel")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 25,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ EDU_REC_2 * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_edu[prev_year_edu$subgroup=="Overall",]))

Figure5B<-ggplot(prev_year_edu[prev_year_edu$subgroup=="Subgroups" & prev_year_edu$group=="MARD",], 
       aes(x=YEAR_2, y=prop, group=EDU_REC_2, colour=EDU_REC_2)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Education \nLevel")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 15))+
  geom_text(aes(x = 2,y = 13,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ EDU_REC_2 * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_edu[prev_year_edu$subgroup=="Subgroups" &prev_year_edu$group=="MARD",]))

Figure5C<-ggplot(prev_year_edu[prev_year_edu$subgroup=="Subgroups" & prev_year_edu$group=="MOD",], 
       aes(x=YEAR_2, y=prop, group=EDU_REC_2, colour=EDU_REC_2)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Education \nLevel")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 15))+
  geom_text(aes(x = 2,y = 13,label = 'p value < 0.01'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ EDU_REC_2 * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_edu[prev_year_edu$subgroup=="Subgroups" &prev_year_edu$group=="MOD",]))

Figure5D<-ggplot(prev_year_edu[prev_year_edu$subgroup=="Subgroups" & prev_year_edu$group=="SIRD",], 
       aes(x=YEAR_2, y=prop, group=EDU_REC_2, colour=EDU_REC_2)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Education \nLevel")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 10,label = 'p value = <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ EDU_REC_2 * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_edu[prev_year_edu$subgroup=="Subgroups" &prev_year_edu$group=="SIRD",]))

Figure5E<-ggplot(prev_year_edu[prev_year_edu$subgroup=="Subgroups" & prev_year_edu$group=="SIID",], 
       aes(x=YEAR_2, y=prop, group=EDU_REC_2, colour=EDU_REC_2)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Education \nLevel")+ggtitle("SIDD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 10,label = 'p value = <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ EDU_REC_2 * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_edu[prev_year_edu$subgroup=="Subgroups" &prev_year_edu$group=="SIID",]))

Figure5_Down<-ggpubr::ggarrange(Figure5B,Figure5C,Figure5D,Figure5E,nrow = 2, ncol = 2,common.legend = T)
Figure5<-ggpubr::ggarrange(Figure5A,Figure5_Down,nrow = 2, ncol = 1,labels = c("C"),common.legend = T)

Figure2_FINAL<-ggpubr::ggarrange(Figure2,Figure3,Figure5,nrow = 1, ncol = 3,common.legend = T)

ggsave(Figure2_FINAL,
       filename = "Figure2.jpg", 
       width = 50, 
       height = 28,
       units=c("cm"),
       dpi = 400,
       limitsize = FALSE)

#----Sex and etnicity stratification-----

DM2_prev_race_sex<-svyby(~DM2_final,~YEAR_2+RACE_REC+GENDER_REC, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_race_sex<-svyby(~ClusterMARD, ~YEAR_2+RACE_REC+GENDER_REC, NHANES_subset, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_race_sex<-svyby(~ClusterMOD, ~YEAR_2+RACE_REC+GENDER_REC , NHANES_subset, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_race_sex<-svyby(~ClusterSIRD,~YEAR_2+RACE_REC+GENDER_REC, NHANES_subset, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_race_sex<-svyby(~ClusterSIID, ~YEAR_2+RACE_REC+GENDER_REC, NHANES_subset, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_race_sex)[names(DM2_prev_race_sex) == "DM2_final"] <- "prop"
names(DM2_MARD_race_sex)[names(DM2_MARD_race_sex) == "ClusterMARD"] <- "prop"
names(DM2_MOD_race_sex)[names(DM2_MOD_race_sex) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_race_sex)[names(DM2_SIRD_race_sex) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_race_sex)[names(DM2_SIID_race_sex) == "ClusterSIID"] <- "prop"

DM2_prev_race_sex$group<-"Type 2 Diabetes"
DM2_MARD_race_sex$group<-"MARD"
DM2_MOD_race_sex$group<-"MOD"
DM2_SIRD_race_sex$group<-"SIRD"
DM2_SIID_race_sex$group<-"SIID"

prev_year_race_sex<-rbind(DM2_prev_race_sex,DM2_MARD_race_sex,DM2_MOD_race_sex,DM2_SIRD_race_sex,DM2_SIID_race_sex)
prev_year_race_sex$subgroup<-ifelse(prev_year_race_sex$group=="Type 2 Diabetes", 0, 1)
prev_year_race_sex$subgroup<-factor(prev_year_race_sex$subgroup, labels = c("Overall","Subgroups"))

Supp_Fig_X<-ggplot(prev_year_race_sex, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~RACE_REC+GENDER_REC+subgroup, nrow=3, scales = "free")

ggsave(Supp_Fig_X,
       filename = "Supp_Fig_11.jpg", 
       width = 75, 
       height = 20,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)


#----Obesity weighted prevalence-----

NHANES_DM2$OBESITY_2<-NULL;
NHANES_DM2$OBESITY_2[NHANES_DM2$BMI_2<30]<-0
NHANES_DM2$OBESITY_2[NHANES_DM2$BMI_2>=30]<-1

NHANES_DM2_ob<-NHANES_DM2%>%filter(BMI_2>=0)

NHANES_all_ob <- svydesign(data=NHANES_DM2_ob, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
NHANES_subset_ob <- subset(NHANES_all_ob, BMI_2>=0)
Ob_prev<-svyby(~OBESITY_2, ~YEAR_2, NHANES_subset_ob, svymean) %>%   mutate(OBESITY_2 = round(OBESITY_2*100, digits=1), se=round(se*100, digits=1))

#----Ethnicity and obesity stratification-----

DM2_prev_race_obesity<-svyby(~DM2_final,~YEAR_2+RACE_REC+BMI_CAT, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_race_obesity<-svyby(~ClusterMARD, ~YEAR_2+RACE_REC+BMI_CAT, NHANES_subset, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_race_obesity<-svyby(~ClusterMOD, ~YEAR_2+RACE_REC+BMI_CAT , NHANES_subset, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_race_obesity<-svyby(~ClusterSIRD,~YEAR_2+RACE_REC+BMI_CAT, NHANES_subset, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_race_obesity<-svyby(~ClusterSIID, ~YEAR_2+RACE_REC+BMI_CAT, NHANES_subset, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_race_obesity)[names(DM2_prev_race_obesity) == "DM2_final"] <- "prop"
names(DM2_MARD_race_obesity)[names(DM2_MARD_race_obesity) == "ClusterMARD"] <- "prop"
names(DM2_MOD_race_obesity)[names(DM2_MOD_race_obesity) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_race_obesity)[names(DM2_SIRD_race_obesity) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_race_obesity)[names(DM2_SIID_race_obesity) == "ClusterSIID"] <- "prop"

DM2_prev_race_obesity$group<-"Type 2 Diabetes"
DM2_MARD_race_obesity$group<-"MARD"
DM2_MOD_race_obesity$group<-"MOD"
DM2_SIRD_race_obesity$group<-"SIRD"
DM2_SIID_race_obesity$group<-"SIID"

prev_year_race_obesity<-rbind(DM2_prev_race_obesity,DM2_MARD_race_obesity,DM2_MOD_race_obesity,DM2_SIRD_race_obesity,DM2_SIID_race_obesity)
prev_year_race_obesity$subgroup<-ifelse(prev_year_race_obesity$group=="Type 2 Diabetes", 0, 1)
prev_year_race_obesity$subgroup<-factor(prev_year_race_obesity$subgroup, labels = c("Overall","Subgroups"))


#----Etnicity, age and obesity stratification-----

DM2_prev_race_age_obesity<-svyby(~DM2_final,~YEAR_2+RACE_REC+age_cat+BMI_CAT, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_race_age_obesity<-svyby(~ClusterMARD, ~YEAR_2+RACE_REC+age_cat+BMI_CAT, NHANES_subset, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_race_age_obesity<-svyby(~ClusterMOD, ~YEAR_2+RACE_REC+age_cat+BMI_CAT , NHANES_subset, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_race_age_obesity<-svyby(~ClusterSIRD,~YEAR_2+RACE_REC+age_cat+BMI_CAT, NHANES_subset, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_race_age_obesity<-svyby(~ClusterSIID, ~YEAR_2+RACE_REC+age_cat+BMI_CAT, NHANES_subset, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_race_age_obesity)[names(DM2_prev_race_age_obesity) == "DM2_final"] <- "prop"
names(DM2_MARD_race_age_obesity)[names(DM2_MARD_race_age_obesity) == "ClusterMARD"] <- "prop"
names(DM2_MOD_race_age_obesity)[names(DM2_MOD_race_age_obesity) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_race_age_obesity)[names(DM2_SIRD_race_age_obesity) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_race_age_obesity)[names(DM2_SIID_race_age_obesity) == "ClusterSIID"] <- "prop"

DM2_prev_race_age_obesity$group<-"Type 2 Diabetes"
DM2_MARD_race_age_obesity$group<-"MARD"
DM2_MOD_race_age_obesity$group<-"MOD"
DM2_SIRD_race_age_obesity$group<-"SIRD"
DM2_SIID_race_age_obesity$group<-"SIID"

prev_year_race_age_obesity<-rbind(DM2_prev_race_age_obesity,DM2_MARD_race_age_obesity,DM2_MOD_race_age_obesity,DM2_SIRD_race_age_obesity,DM2_SIID_race_age_obesity)
prev_year_race_age_obesity$subgroup<-ifelse(prev_year_race_age_obesity$group=="Type 2 Diabetes", 0, 1)
prev_year_race_age_obesity$subgroup<-factor(prev_year_race_age_obesity$subgroup, labels = c("Overall","Subgroups"))


#----Trend analysis-----

#Overall diabetes and subgroups trends

DM2_prev[,2]
DM2_prev[,2]-DM2_prev[,3]
DM2_prev[,2]+DM2_prev[,3]
summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_prev))

DM2_MARD[,2]
DM2_MARD[,2]-DM2_MARD[,3]
DM2_MARD[,2]+DM2_MARD[,3]
summary(glm(prop~as.numeric(YEAR_2), data=DM2_MARD))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_MARD))

DM2_MOD[,2]
DM2_MOD[,2]-DM2_MOD[,3]
DM2_MOD[,2]+DM2_MOD[,3]
summary(glm(prop~as.numeric(YEAR_2), data=DM2_MOD))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_MOD))

DM2_SIRD[,2]
DM2_SIRD[,2]-DM2_SIRD[,3]
DM2_SIRD[,2]+DM2_SIRD[,3]
summary(glm(prop~as.numeric(YEAR_2), data=DM2_SIRD))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_SIRD))

DM2_SIID[,2]
DM2_SIID[,2]-DM2_SIID[,3]
DM2_SIID[,2]+DM2_SIID[,3]
summary(glm(prop~as.numeric(YEAR_2), data=DM2_SIID))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_SIID))

#Overall diabetes trends in sex

#Overall diabetes trend by etnicities

View(DM2_prev_sex)
head(DM2_prev_sex)

summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev_sex, subset = GENDER_REC=="Men"))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_prev_sex, subset = GENDER_REC=="Men"))

summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev_sex, subset = GENDER_REC=="Women"))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_prev_sex, subset = GENDER_REC=="Women"))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men"  & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men"  & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & group=="SIID")))


#Subgroups by etnicities

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men" & RACE_REC =="Non-Hispanic Black" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men" & RACE_REC =="Non-Hispanic Black" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic Black" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic Black" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic Black" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic Black" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic Black" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic Black" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men" & RACE_REC =="Non-Hispanic White" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men" & RACE_REC =="Non-Hispanic White" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic White" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic White" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic White" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic White" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic White" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Non-Hispanic White" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men" & RACE_REC =="Mexican American" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Men" & RACE_REC =="Mexican American" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Mexican American" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Mexican American" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Mexican American" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Mexican American" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Mexican American" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Men" & RACE_REC =="Mexican American" & group=="SIID")))




summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Women" & RACE_REC =="Non-Hispanic Black" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Women" & RACE_REC =="Non-Hispanic Black" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic Black" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic Black" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic Black" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic Black" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic Black" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic Black" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Women" & RACE_REC =="Non-Hispanic White" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Women" & RACE_REC =="Non-Hispanic White" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic White" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic White" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic White" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic White" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic White" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Non-Hispanic White" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Women" & RACE_REC =="Mexican American" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC =="Women" & RACE_REC =="Mexican American" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Mexican American" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Mexican American" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Mexican American" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Mexican American" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Mexican American" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_sex, subset = (subgroup=="Subgroups" & GENDER_REC=="Women" & RACE_REC =="Mexican American" & group=="SIID")))


#Overall diabetes trend by etnicities

summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev_race, subset = RACE_REC=="Non-Hispanic Black"))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_prev_race, subset = RACE_REC=="Non-Hispanic Black"))

summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev_race, subset = RACE_REC=="Non-Hispanic White"))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_prev_race, subset = RACE_REC=="Non-Hispanic White"))

summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev_race, subset = RACE_REC=="Mexican American"))
confint(glm(prop~as.numeric(YEAR_2), data=DM2_prev_race, subset = RACE_REC=="Mexican American"))

#Subgroups by etnicities
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & group=="SIID")))


#Overall diabetes trend by education 
table(prev_year_edu$EDU_REC_2)

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Primary School \nor \nNo-Education")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Primary School \nor \nNo-Education")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="College \nor \nHigher")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="College \nor \nHigher")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="College \nor \nHigher" & group=="SIID")))


#Education and etinicity stratification
table(prev_year_race_edu$RACE_REC)
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="College \nor \nHigher" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="College \nor \nHigher" & group=="SIRD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="College \nor \nHigher" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="College \nor \nHigher" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Primary School \nor \nNo-Education" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="College \nor \nHigher" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="College \nor \nHigher" & group=="SIID")))





#Diabetes age stratification

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Overall" & age_cat=="≥18-39 years")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Overall" & age_cat=="≥18-39 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Overall" & age_cat=="40-64 years")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Overall" & age_cat=="40-64 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Overall" & age_cat=="≥65 years")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Overall" & age_cat=="≥65 years")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="MARD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="MOD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="SIRD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥18-39 years" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="40-64 years" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_age, subset = (subgroup=="Subgroups" & age_cat=="≥65 years" & group=="SIID")))

#Diabetes BMI stratification

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Normal-Weight")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Normal-Weight")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Overweight")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Overweight")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Obese")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Obese")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="MARD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="MOD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="SIRD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Normal-Weight" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Overweight" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & BMI_CAT=="Obese" & group=="SIID")))


#Diabetes age and etnicity stratification

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="≥18-39 years" & RACE_REC=="Non-Hispanic Black")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="40-64 years" & RACE_REC=="Non-Hispanic Black")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="≥65 years" & RACE_REC=="Non-Hispanic Black")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="≥18-39 years" & RACE_REC=="Mexican American")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="40-64 years" & RACE_REC=="Mexican American")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="≥65 years" & RACE_REC=="Mexican American")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="≥18-39 years" & RACE_REC=="Non-Hispanic White")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="40-64 years" & RACE_REC=="Non-Hispanic White")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Overall" & age_cat=="≥65 years" & RACE_REC=="Non-Hispanic White")))

#Trend in subgroups for age categories

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & age_cat=="≥18-39 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & age_cat=="40-64 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & age_cat=="≥65 years")))

#Trend in subgroups for age and etnicity categories

#18-39 years 
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Mexican American" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Non-Hispanic White" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥18-39 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Mexican American" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Non-Hispanic White" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥18-39 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Mexican American" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Non-Hispanic White" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥18-39 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Mexican American" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Non-Hispanic White" & age_cat=="≥18-39 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥18-39 years")))

#40-64 years
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Mexican American" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Non-Hispanic White" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Non-Hispanic Black" & age_cat=="40-64 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Mexican American" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Non-Hispanic White" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Non-Hispanic Black" & age_cat=="40-64 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Mexican American" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Non-Hispanic White" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Non-Hispanic Black" & age_cat=="40-64 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Mexican American" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Non-Hispanic White" & age_cat=="40-64 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Non-Hispanic Black" & age_cat=="40-64 years")))

#≥65 years
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Mexican American" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Non-Hispanic White" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MARD" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥65 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Mexican American" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Non-Hispanic White" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="MOD" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥65 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Mexican American" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Non-Hispanic White" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIRD" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥65 years")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Mexican American" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Non-Hispanic White" & age_cat=="≥65 years")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_age, subset = (subgroup=="Subgroups" & group=="SIID" & RACE_REC=="Non-Hispanic Black" & age_cat=="≥65 years")))

#Stratification for BMI cathegories 

View(prev_year_obesity)
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Normal-Weight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Overweight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Overall" & BMI_CAT=="Obese")))


#Trend in subgroups for age categories


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="MARD" & BMI_CAT=="Normal-Weight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="MOD" & BMI_CAT=="Normal-Weight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="SIRD" & BMI_CAT=="Normal-Weight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="SIID" & BMI_CAT=="Normal-Weight")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="MARD" & BMI_CAT=="Overweight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="MOD" & BMI_CAT=="Overweight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="SIRD" & BMI_CAT=="Overweight")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="SIID" & BMI_CAT=="Overweight")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="MARD" & BMI_CAT=="Obese")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="MOD" & BMI_CAT=="Obese")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="SIRD" & BMI_CAT=="Obese")))
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_obesity, subset = (subgroup=="Subgroups" & group=="SIID" & BMI_CAT=="Obese")))

#----SAID stratification in NHANES III (Sup Fig 2)-----

NHANES_DM2_III<-NHANES_DM2%>%dplyr::select(seqn,Cluster,YEAR)%>%na.omit()

nhanes_3<-read_csv("nhanes_3_diabetes.csv")
nhanes_3$Glucose[nhanes_3$Glucose == 88888] <- NA
nhanes_3$Glycohemoglobin[nhanes_3$Glycohemoglobin == 8888] <- NA
nhanes_3$DM2_final<-as.numeric((nhanes_3$Glycohemoglobin>=6.5 |nhanes_3$Glucose>=126))
nhanes_3$DM2_final<-na.tools::na.replace(nhanes_3$DM2_final,0)

summary(nhanes_3$GAD_65)

nhanes_3$seqn<-paste0(str_pad(nhanes_3$seqn, 1,pad = "0"),c("-3"))
nhanes_3<-nhanes_3%>%left_join(NHANES_DM2_III,by="seqn")


nhanes_3$Cluster_REC<-NULL
nhanes_3$Cluster_REC[nhanes_3$Cluster=="MARD"]<-"MARD"
nhanes_3$Cluster_REC[nhanes_3$Cluster=="MOD"]<-"MOD"
nhanes_3$Cluster_REC[nhanes_3$Cluster=="SIID"]<-"SIID"
nhanes_3$Cluster_REC[nhanes_3$Cluster=="SIRD"]<-"SIRD"
nhanes_3$Cluster_REC[nhanes_3$DM2_final==1 & nhanes_3$GAD_65>=0.069]<-"SAID"

nhanes_3$AGE_18<-NULL;nhanes_3$AGE_18[nhanes_3$AGE>=18]<-1;nhanes_3$AGE_18[nhanes_3$AGE<18]<-0
nhanes_3$BMI[nhanes_3$BMI == 8888] <- NA
nhanes_3$Weight[nhanes_3$Weight == 888888] <- NA
nhanes_3$Height[nhanes_3$Height == 88888.0] <- NA
nhanes_3$BMI_2<-nhanes_3$Weight/((nhanes_3$Height/100)^2)
nhanes_3<-nhanes_3%>%mutate(one = 1)

d1<-dummies::dummy(nhanes_3$Cluster_REC)
nhanes_3<-cbind(nhanes_3,d1)

nhanes_3$RACE_REC<-NULL
nhanes_3$RACE_REC[nhanes_3$RACE==1]<-2
nhanes_3$RACE_REC[nhanes_3$RACE==2]<-1
nhanes_3$RACE_REC[nhanes_3$RACE==3]<-3
nhanes_3$RACE_REC[nhanes_3$RACE==4]<-2
nhanes_3$RACE_REC<-factor(nhanes_3$RACE_REC, labels = c("Mexican American","Non-Hispanic White", "Non-Hispanic Black"))

nhanes_3$age_cat<-NULL;
nhanes_3$age_cat[nhanes_3$AGE>=0 & nhanes_3$AGE<=17]<-0
nhanes_3$age_cat[nhanes_3$AGE>=18 & nhanes_3$AGE<=39]<-1
nhanes_3$age_cat[nhanes_3$AGE>=40 & nhanes_3$AGE<=64]<-2
nhanes_3$age_cat[nhanes_3$AGE>=65]<-3
nhanes_3$age_cat<-factor(nhanes_3$age_cat,labels = c("<18 years", "≥18-39 years", "40-64 years","≥65 years"))

nhanes_3$OBESITY_2<-NULL;
nhanes_3$OBESITY_2[nhanes_3$BMI_2<30]<-0
nhanes_3$OBESITY_2[nhanes_3$BMI_2>=30]<-1
nhanes_3$OBESITY_2<-factor(nhanes_3$OBESITY_2, labels = c("Non-Obese", "Obese"))

nhanes_3$BMI_CAT<-NULL;
nhanes_3$BMI_CAT[nhanes_3$BMI_2<25]<-0
nhanes_3$BMI_CAT[nhanes_3$BMI_2>=25 & nhanes_3$BMI_2<30]<-1
nhanes_3$BMI_CAT[nhanes_3$BMI_2>=30]<-2
nhanes_3$BMI_CAT<-factor(nhanes_3$BMI_CAT, labels = c("Normal-Weight", "Overweight", "Obese"))

nhanes_3$Cluster_RECMARD
nhanes_3$Cluster_RECMOD
nhanes_3$Cluster_RECSIRD
nhanes_3$Cluster_RECSIID
nhanes_3$Cluster_RECSAID

NHANES_all_3 <- svydesign(data=nhanes_3, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
NHANES_subset_3 <- subset(NHANES_all_3, AGE_18==1)
NHANES_subset1_3<- subset(NHANES_all_3, AGE_18==1 & !is.na(nhanes_3$Insulin))

DM2_prev<-svyby(~DM2_final, ~YEAR.x,NHANES_subset_3, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD<-svyby(~Cluster_RECMARD, ~YEAR.x,NHANES_subset1_3, svymean) %>%   mutate(Cluster_RECMARD = round(Cluster_RECMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD<-svyby(~Cluster_RECMOD, ~YEAR.x,NHANES_subset1_3, svymean) %>%   mutate(Cluster_RECMOD = round(Cluster_RECMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD<-svyby(~Cluster_RECSIRD, ~YEAR.x,NHANES_subset1_3, svymean) %>%  mutate(Cluster_RECSIRD = round(Cluster_RECSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID<-svyby(~Cluster_RECSIID, ~YEAR.x,NHANES_subset1_3, svymean) %>%  mutate(Cluster_RECSIID = round(Cluster_RECSIID*100, digits=1), se=round(se*100, digits=1))
DM2_SAID<-svyby(~Cluster_RECSAID, ~YEAR.x,NHANES_subset1_3, svymean) %>%  mutate(Cluster_RECSAID = round(Cluster_RECSAID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev)[names(DM2_prev) == "DM2_final"] <- "prop"
names(DM2_MARD)[names(DM2_MARD) == "Cluster_RECMARD"] <- "prop"
names(DM2_MOD)[names(DM2_MOD) == "Cluster_RECMOD"] <- "prop"
names(DM2_SIRD)[names(DM2_SIRD) == "Cluster_RECSIRD"] <- "prop"
names(DM2_SIID)[names(DM2_SIID) == "Cluster_RECSIID"] <- "prop"
names(DM2_SAID)[names(DM2_SAID) == "Cluster_RECSAID"] <- "prop"

DM2_prev$group<-"Type 2 Diabetes"
DM2_MARD$group<-"MARD"
DM2_MOD$group<-"MOD"
DM2_SIRD$group<-"SIRD"
DM2_SIID$group<-"SIDD"
DM2_SAID$group<-"SAID"

prev_year<-rbind(DM2_prev,DM2_MARD,DM2_MOD,DM2_SIRD,DM2_SIID,DM2_SAID)
prev_year$subgroup<-ifelse(prev_year$group=="Type 2 Diabetes", 1, 0)
prev_year$subgroup<-factor(prev_year$subgroup, labels = c("Subgroups", "Overall"))


Supp_Fig_6A<-ggplot(prev_year, aes(x=YEAR.x, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  scale_color_jama()+
  labs(colour="Diabetes Clusters")+facet_wrap(~subgroup, scales = "free")+
  scale_y_continuous(limits = c(0,8))+
  theme(legend.position="none")


NHANES_DM2_Prev_III<-nhanes_3%>%dplyr::select(Cluster_REC,YEAR.x)%>%na.omit()

NHANES_DM2_Prev_FIG<-NHANES_DM2_Prev_III%>%
  group_by(YEAR.x,Cluster_REC)%>%
  summarise(freq_clus=n())%>%
  mutate(prop_clus_year=(freq_clus/sum(freq_clus))*100)%>%
  na.omit()

NHANES_DM2_Prev_FIG$Cluster_REC<-factor(NHANES_DM2_Prev_FIG$Cluster_REC,labels = c("MARD","MOD","SAID","SIDD","SIRD"))

Supp_Fig_6B<-ggplot(NHANES_DM2_Prev_FIG, aes(y=prop_clus_year, x=Cluster_REC,fill=Cluster_REC)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(width=0.8))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme_hc()+
  ylab("Absolute Proportion (%)")+
  geom_text(position = position_dodge(width= 0.8),aes(label=round(prop_clus_year*1,2)), vjust=1.6,color="white", size=3)+
  xlab("NHANES Cycle")+
  scale_fill_jama()+
  labs(fill="Diabetes Clusters")+
  facet_wrap(~YEAR.x, nrow = 1)

NHANES_DM2_Cluster_III<-nhanes_3%>%dplyr::select(Cluster_REC,Glycohemoglobin,HOMA2_IR,HOMA2_B,Age_diabetes,BMI,BMI_CAT,RACE_REC)%>%na.omit()
my_comparisons <- list( c("MARD", "MOD"), c("MARD", "SIID"), c("MARD", "SIRD"), c("MARD", "SAID"), 
                        c("MOD", "SIID"), c("MOD", "SIRD"), c("MOD", "SAID"),
                        c("SIID", "SIRD"), c("SIID", "SAID"))

NHANES_DM2_Cluster_III$Age_diabetes[NHANES_DM2_Cluster_III$Age_diabetes == 999] <- NA
NHANES_DM2_Cluster_III$Age_diabetes[NHANES_DM2_Cluster_III$Age_diabetes == 888] <- NA
NHANES_DM2_Cluster_III$Cluster_REC<-factor(NHANES_DM2_Cluster_III$Cluster_REC,labels = c("MARD","MOD","SAID","SIDD","SIRD"))
#Frecuncy of clusters
data<-as.data.frame(prop.table(table(NHANES_DM2_Cluster_III$Cluster_REC))*100)
names(data)<-c("Diabetes Cluster", "Absolute Frecuency")

data <- data %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

data$abs_freq_2<-as.numeric(format(round(data$`Absolute Frecuency`, 2), nsmall = 2))
data$`Diabetes Cluster`<-factor(data$`Diabetes Cluster`,labels = c("MARD","MOD","SAID","SIDD","SIRD"))

Supp_Fig_6C<-ggplot(data, aes(x="", y=`Absolute Frecuency`, fill=`Diabetes Cluster`)) +
  geom_bar(stat="identity", width=2) +
  coord_polar("y", start=0)+
  theme_void()+
  ylab("")+
  theme(legend.position="bottom") +
  geom_text(aes(y = ypos, label = abs_freq_2), color = "white", size=6) +
  scale_fill_jama()

#Hba1c

Supp_Fig_6D<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=Glycohemoglobin, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Glycated hemoglobin (%)")+
  scale_fill_jama()+
  theme(legend.position="none")

#HOMA2-IR
Sup_Fig_6E<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=HOMA2_IR, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-IR")+
  scale_fill_jama()+
  scale_y_log10()+  
  theme(legend.position="none")

#HOMA2-B
Sup_Fig_6F<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=HOMA2_B, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-B")+
  scale_fill_jama()+
  scale_y_log10()+  
  theme(legend.position="none")

#Age at diagnosis
Sup_Fig_6G<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=Age_diabetes, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Age at Diabetes Diagnosis
(years)")+
  scale_fill_jama()+
  theme(legend.position="none")

#Body Mass Index
Sup_Fig_6H<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=BMI, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Body Mass Index (kg/m2)")+
  scale_fill_jama()+
  theme(legend.position="none")


Sup_Fig_6_Upp_R<-ggarrange(Supp_Fig_6D,Sup_Fig_6E,Sup_Fig_6F,labels= c("C","D","E"),ncol = 3, nrow = 1)
Sup_Fig_6_Down_R<-ggarrange(Sup_Fig_6G,Sup_Fig_6H,labels= c("F","G"),ncol = 2, nrow = 1)
Sup_Fig_6_LOW_RIGHT<-ggarrange(Sup_Fig_6_Upp_R,Sup_Fig_6_Down_R, nrow = 2,ncol = 1)
Sup_Fig_6_DOWN<-ggarrange(Sup_Fig_6_LOW_RIGHT,nrow = 1,ncol = 1,common.legend = T)
Sup_Fig_6_UP<-ggarrange(Supp_Fig_6A,Supp_Fig_6B, labels= c("A","B"),nrow = 1,ncol = 2,common.legend = T)
Sup_Fig_8<-ggarrange(Sup_Fig_6_UP,Sup_Fig_6_DOWN, nrow = 2,ncol = 1,common.legend = T)

ggsave(Sup_Fig_8,
       filename = "Supp_Fig_2.jpg", 
       width = 45, 
       height = 30,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)


#----Differences of diabetes (Sup Fig 3 & 4)-----

NHANES_DM2_Cluster<-NHANES_DM2%>%dplyr::select(Cluster,Glycohemoglobin,HOMA2_IR,HOMA2_B,Age_diabetes,BMI,BMI_CAT,RACE_REC)%>%na.omit()

my_comparisons <- list( c("MARD", "MOD"), c("MARD", "SIID"), c("MARD", "SIRD"), 
                        c("MOD", "SIID"), c("MOD", "SIRD"),
                        c("SIID", "SIRD"))

#Frecuncy of clusters
data<-as.data.frame(prop.table(table(NHANES_DM2_Cluster$Cluster))*100)
names(data)<-c("Diabetes Cluster", "Absolute Frecuency")

data <- data %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

data$abs_freq_2<-as.numeric(format(round(data$`Absolute Frecuency`, 2), nsmall = 2))
data$`Diabetes Cluster`<-factor(data$`Diabetes Cluster`,labels = c("MARD","MOD","SIDD","SIRD"))

Sup_Fig_3A<-ggplot(data, aes(x="", y=`Absolute Frecuency`, fill=`Diabetes Cluster`)) +
  geom_bar(stat="identity", width=2) +
  coord_polar("y", start=0)+
  theme_void()+
  ylab("")+
  theme(legend.position="bottom") +
  geom_text(aes(y = ypos, label = abs_freq_2), color = "white", size=6) +
  scale_fill_jama()

#Hba1c
NHANES_DM2_Cluster$Cluster<-factor(NHANES_DM2_Cluster$Cluster,labels = c("MARD","MOD","SIDD","SIRD"))
Sup_Fig_3B<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=Glycohemoglobin, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Glycated hemoglobin (%)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#HOMA2-IR
Sup_Fig_3C<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=HOMA2_IR, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-IR")+
  scale_fill_jama()+
  scale_y_log10()+  
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#HOMA2-B
Sup_Fig_3D<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=HOMA2_B, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-B")+
  scale_fill_jama()+
  scale_y_log10()+  
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#Age at diagnosis
Sup_Fig_3E<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=Age_diabetes, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Age at Diabetes Diagnosis
(years)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#Body Mass Index
Sup_Fig_3F<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=BMI, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Body Mass Index (kg/m2)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

Sup_Fig_3_Upp_R<-ggarrange(Sup_Fig_3B,Sup_Fig_3C,Sup_Fig_3D,ncol = 3, nrow = 1, labels = c("B","C","D"))
Sup_Fig_3_Down_R<-ggarrange(Sup_Fig_3E,Sup_Fig_3F,ncol = 2, nrow = 1, labels = c("E","F"))
Sup_Fig_3_R<-ggarrange(Sup_Fig_3_Upp_R,Sup_Fig_3_Down_R, nrow = 2,ncol = 1)
Sup_Fig_6<-ggarrange(Sup_Fig_3A,Sup_Fig_3_R, nrow = 1,ncol = 2, labels = c("A"), common.legend = T)

ggsave(Sup_Fig_6,
       filename = "Supp_Fig_3.jpg", 
       width = 40, 
       height = 20,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

#Etnicity stratification

data1<-as.data.frame(prop.table(table(NHANES_DM2_Cluster[NHANES_DM2_Cluster$RACE_REC=="Mexican American",]$Cluster))*100)
names(data1)<-c("Diabetes Cluster", "Absolute Frecuency")

data1 <- data1 %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
data1$abs_freq_2<-as.numeric(format(round(data1$`Absolute Frecuency`, 2), nsmall = 2))
data1$`Diabetes Cluster`<-factor(data1$`Diabetes Cluster`,labels = c("MARD","MOD","SIDD","SIRD"))


data2<-as.data.frame(prop.table(table(NHANES_DM2_Cluster[NHANES_DM2_Cluster$RACE_REC=="Non-Hispanic White",]$Cluster))*100)
names(data2)<-c("Diabetes Cluster", "Absolute Frecuency")

data2 <- data2 %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
data2$abs_freq_2<-as.numeric(format(round(data2$`Absolute Frecuency`, 2), nsmall = 2))
data2$`Diabetes Cluster`<-factor(data2$`Diabetes Cluster`,labels = c("MARD","MOD","SIDD","SIRD"))

data3<-as.data.frame(prop.table(table(NHANES_DM2_Cluster[NHANES_DM2_Cluster$RACE_REC=="Non-Hispanic Black",]$Cluster))*100)
names(data3)<-c("Diabetes Cluster", "Absolute Frecuency")

data3 <- data3 %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
data3$abs_freq_2<-as.numeric(format(round(data3$`Absolute Frecuency`, 2), nsmall = 2))
data3$`Diabetes Cluster`<-factor(data3$`Diabetes Cluster`,labels = c("MARD","MOD","SIDD","SIRD"))

Sup_Fig_4A.1<-ggplot(data1, aes(x="", y=`Absolute Frecuency`, fill=`Diabetes Cluster`)) +
  geom_bar(stat="identity", width=2) +
  coord_polar("y", start=0)+
  theme_void()+
  ylab("")+
  theme(legend.position="bottom") +
  geom_text(aes(y = ypos, label = abs_freq_2), color = "white", size=6) +
  scale_fill_jama()+
  ggtitle("Mexican American")

Sup_Fig_4A.2<-ggplot(data2, aes(x="", y=`Absolute Frecuency`, fill=`Diabetes Cluster`)) +
  geom_bar(stat="identity", width=2) +
  coord_polar("y", start=0)+
  theme_void()+
  ylab("")+
  theme(legend.position="bottom") +
  geom_text(aes(y = ypos, label = abs_freq_2), color = "white", size=6) +
  scale_fill_jama()+
  ggtitle("Non-Hispanic White")

Sup_Fig_4A.3<-ggplot(data3, aes(x="", y=`Absolute Frecuency`, fill=`Diabetes Cluster`)) +
  geom_bar(stat="identity", width=2) +
  coord_polar("y", start=0)+
  theme_void()+
  ylab("")+
  theme(legend.position="bottom") +
  geom_text(aes(y = ypos, label = abs_freq_2), color = "white", size=6) +
  scale_fill_jama()+
  ggtitle("Non-Hispanic Black")

Sup_Fig_4_L<-ggarrange(Sup_Fig_4A.1,Sup_Fig_4A.2,Sup_Fig_4A.3,ncol = 3, nrow = 1, labels = c("A"), common.legend = T)

#Hba1c

Sup_Fig_4B<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=Glycohemoglobin, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Glycated hemoglobin (%)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")+
  facet_wrap(~RACE_REC)

#HOMA2-IR
Sup_Fig_4C<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=HOMA2_IR, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-IR")+
  scale_fill_jama()+
  scale_y_log10()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")+
  facet_wrap(~RACE_REC)

#HOMA2-B
Sup_Fig_4D<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=HOMA2_B, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-B")+
  scale_fill_jama()+
  scale_y_log10()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")+
  facet_wrap(~RACE_REC)

#Age at diagnosis

Sup_Fig_4E<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=Age_diabetes, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Age at Diabetes Diagnosis(years)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")+
  facet_wrap(~RACE_REC)

#Body Mass Index
Sup_Fig_4F<-ggplot(NHANES_DM2_Cluster, aes(x=Cluster,y=BMI, fill=Cluster))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Body Mass Index (kg/m2)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")+
  facet_wrap(~RACE_REC)

Sup_Fig_4_Upp_R<-ggarrange(Sup_Fig_4B,Sup_Fig_4C,Sup_Fig_4D,ncol = 3, nrow = 1, labels = c("B","C","D"))
Sup_Fig_4_Down_R<-ggarrange(Sup_Fig_4E,Sup_Fig_4F,ncol = 2, nrow = 1, labels = c("E","F"))
Sup_Fig_4_R<-ggarrange(Sup_Fig_4_Upp_R,Sup_Fig_4_Down_R, nrow = 2,ncol = 1)
Sup_Fig_7<-ggarrange(Sup_Fig_4_L,Sup_Fig_4_R, nrow = 1,ncol = 2, common.legend = T)

ggsave(Sup_Fig_7,
       filename = "Supp_Fig_4.jpg", 
       width = 80, 
       height = 18,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

#----Proportion of diabetes in the studied period (Sup Fig 5)-----
NHANES_DM2$YEAR_2
NHANES_DM2_Prev<-NHANES_DM2%>%filter(!is.na(Cluster))%>%
  dplyr::select(Cluster,YEAR_2)%>%na.omit()

NHANES_DM2_Prev_FIG<-NHANES_DM2_Prev%>%
  group_by(YEAR_2,Cluster)%>%
  summarise(freq_clus=n())%>%
  mutate(prop_clus_year=(freq_clus/sum(freq_clus))*100)%>%
  na.omit()

Supp_Fig_5<-ggplot(NHANES_DM2_Prev_FIG, aes(x=YEAR_2, y=prop_clus_year, group=Cluster, colour=Cluster)) + 
  geom_line(size=2.0) +
  geom_point(size=2.5)+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Diabetes Subgroups")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_continuous(limits = c(10,40))

ggsave(Supp_Fig_5,
       filename = "Supp_Fig_5.jpg", 
       width = 20, 
       height = 15,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

#----Education and etnicity stratification (Sup Fig 6)----

DM2_prev_race_edu<-svyby(~DM2_final,~YEAR_2+RACE_REC+EDU_REC_2, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_race_edu<-svyby(~ClusterMARD, ~YEAR_2+RACE_REC+EDU_REC_2, NHANES_subset, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_race_edu<-svyby(~ClusterMOD, ~YEAR_2+RACE_REC+EDU_REC_2 , NHANES_subset, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_race_edu<-svyby(~ClusterSIRD,~YEAR_2+RACE_REC+EDU_REC_2, NHANES_subset, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_race_edu<-svyby(~ClusterSIID, ~YEAR_2+RACE_REC+EDU_REC_2, NHANES_subset, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_race_edu)[names(DM2_prev_race_edu) == "DM2_final"] <- "prop"
names(DM2_MARD_race_edu)[names(DM2_MARD_race_edu) == "ClusterMARD"] <- "prop"
names(DM2_MOD_race_edu)[names(DM2_MOD_race_edu) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_race_edu)[names(DM2_SIRD_race_edu) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_race_edu)[names(DM2_SIID_race_edu) == "ClusterSIID"] <- "prop"

DM2_prev_race_edu$group<-"Type 2 Diabetes"
DM2_MARD_race_edu$group<-"MARD"
DM2_MOD_race_edu$group<-"MOD"
DM2_SIRD_race_edu$group<-"SIRD"
DM2_SIID_race_edu$group<-"SIDD"

prev_year_race_edu<-rbind(DM2_prev_race_edu,DM2_MARD_race_edu,DM2_MOD_race_edu,DM2_SIRD_race_edu,DM2_SIID_race_edu)
prev_year_race_edu$subgroup<-ifelse(prev_year_race_edu$group=="Type 2 Diabetes", 0, 1)
prev_year_race_edu$subgroup<-factor(prev_year_race_edu$subgroup, labels = c("Overall","Subgroups"))

#Primary or No Education
#TYPE 2 DIABETES
Sup_Figure_7A1<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Overall" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("Primary School or No-Education")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 28,label = 'p value = 0.0003'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Overall" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education",]))

#MARD
Sup_Figure_7A2<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="MARD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.0932'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="MARD",]))

#MOD
Sup_Figure_7A3<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="MOD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.109'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="MOD",]))
#SIID
Sup_Figure_7A4<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="SIID",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIID")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.074'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="SIID",]))

#SIRD
Sup_Figure_7A5<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="SIRD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.005'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Primary School \nor \nNo-Education" & prev_year_race_edu$group=="SIRD",]))

#Secondary of High-School
#TYPE 2 DIABETES

Sup_Figure_7B1<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Overall" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("Secondary, High-School or AA Degree")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 28,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Overall" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree",]))

#MARD
Sup_Figure_7B2<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="MARD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4,label = 'p value = 0.0932'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="MARD",]))

#MOD
Sup_Figure_7B3<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="MOD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.109'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="MOD",]))
#SIID
Sup_Figure_7B4<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="SIID",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIID")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.074'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="SIID",]))

#SIRD
Sup_Figure_7B5<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="SIRD",], 
       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.005'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="Secondary, High-School \nor \nAA Degree" & prev_year_race_edu$group=="SIRD",]))


#Collegue or Higher
#TYPE 2 DIABETES

Sup_Figure_7C1<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Overall" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("College or Higher")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 28,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Overall" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher",]))

#MARD
Sup_Figure_7C2<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="MARD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 8))+
  geom_text(aes(x = 2,y = 6,label = 'p value = 0.068'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="MARD",]))

#MOD
Sup_Figure_7C3<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="MOD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 8))+
  geom_text(aes(x = 2,y = 6,label = 'p value = 0.068'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="MOD",]))
#SIID
Sup_Figure_7C4<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="SIID",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIID")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value = 0.0013'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="SIID",]))

#SIRD
Sup_Figure_7C5<-ggplot(prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="SIRD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value = 0.373'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_edu[prev_year_race_edu$subgroup=="Subgroups" & prev_year_race_edu$EDU_REC_2=="College \nor \nHigher" & prev_year_race_edu$group=="SIRD",]))

Sup_Figure_7A2_LOW<-ggarrange(Sup_Figure_7A2,Sup_Figure_7A3,Sup_Figure_7A4,Sup_Figure_7A5,ncol = 2,nrow = 2,common.legend = T)
Sup_Figure_7A<-ggarrange(Sup_Figure_7A1,Sup_Figure_7A2_LOW,ncol = 1,nrow = 2,common.legend = T)

Sup_Figure_7B2_LOW<-ggarrange(Sup_Figure_7B2,Sup_Figure_7B3,Sup_Figure_7B4,Sup_Figure_7B5,ncol = 2,nrow = 2,common.legend = T)
Sup_Figure_7B<-ggarrange(Sup_Figure_7B1,Sup_Figure_7B2_LOW,ncol = 1,nrow = 2,common.legend = T)

Sup_Figure_7C2_LOW<-ggarrange(Sup_Figure_7C2,Sup_Figure_7C3,Sup_Figure_7C4,Sup_Figure_7C5,ncol = 2,nrow = 2,common.legend = T)
Sup_Figure_7C<-ggarrange(Sup_Figure_7C1,Sup_Figure_7C2_LOW,ncol = 1,nrow = 2,common.legend = T)

Sup_Figure_7<-ggarrange(Sup_Figure_7A,Sup_Figure_7B,Sup_Figure_7C,ncol = 3,nrow = 1,common.legend = T)

ggsave(Sup_Figure_7,
       filename = "Supp_Fig_6.jpg", 
       width = 50, 
       height = 28,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)


#----Ethnicity and age stratification (Sup Fig 7)-----

DM2_prev_race_age<-svyby(~DM2_final,~YEAR_2+RACE_REC+age_cat, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_race_age<-svyby(~ClusterMARD, ~YEAR_2+RACE_REC+age_cat, NHANES_subset, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_race_age<-svyby(~ClusterMOD, ~YEAR_2+RACE_REC+age_cat , NHANES_subset, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_race_age<-svyby(~ClusterSIRD,~YEAR_2+RACE_REC+age_cat, NHANES_subset, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_race_age<-svyby(~ClusterSIID, ~YEAR_2+RACE_REC+age_cat, NHANES_subset, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_race_age)[names(DM2_prev_race_age) == "DM2_final"] <- "prop"
names(DM2_MARD_race_age)[names(DM2_MARD_race_age) == "ClusterMARD"] <- "prop"
names(DM2_MOD_race_age)[names(DM2_MOD_race_age) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_race_age)[names(DM2_SIRD_race_age) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_race_age)[names(DM2_SIID_race_age) == "ClusterSIID"] <- "prop"

DM2_prev_race_age$group<-"Type 2 Diabetes"
DM2_MARD_race_age$group<-"MARD"
DM2_MOD_race_age$group<-"MOD"
DM2_SIRD_race_age$group<-"SIRD"
DM2_SIID_race_age$group<-"SIDD"

prev_year_race_age<-rbind(DM2_prev_race_age,DM2_MARD_race_age,DM2_MOD_race_age,DM2_SIRD_race_age,DM2_SIID_race_age)
prev_year_race_age$subgroup<-ifelse(prev_year_race_age$group=="Type 2 Diabetes", 0, 1)
prev_year_race_age$subgroup<-factor(prev_year_race_age$subgroup, labels = c("Overall","Subgroups"))


#Primary or No Education
#TYPE 2 DIABETES

Sup_Figure_7A1<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Overall" & prev_year_race_age$age_cat=="≥18-39 years",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("≥18-39 years")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 50))+
  geom_text(aes(x = 2,y = 48,label = 'p value = 0.002'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Overall" & prev_year_race_age$age_cat=="≥18-39 years",]))

#MARD

Sup_Figure_7A2<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="MARD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.0932'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="MARD",]))

#MOD
Sup_Figure_7A3<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="MOD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value = 0.0434'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="MOD",]))
#SIID
Sup_Figure_7A4<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="SIDD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIDD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value = 0.0592'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="SIID",]))

#SIRD
Sup_Figure_7A5<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="SIRD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = ,label = 'p value = 0.139'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥18-39 years" & prev_year_race_age$group=="SIRD",]))

#Secondary of High-School
#TYPE 2 DIABETES

Sup_Figure_7B1<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Overall" & prev_year_race_age$age_cat=="40-64 years",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("40-64 years")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 50))+
  geom_text(aes(x = 2,y = 48,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Overall" & prev_year_race_age$age_cat=="40-64 years",]))

#MARD
Sup_Figure_7B2<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="MARD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.444'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="MARD",]))

#MOD
Sup_Figure_7B3<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="MOD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="MOD",]))
#SIID
Sup_Figure_7B4<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="SIID",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIID")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="SIID",]))

#SIRD
Sup_Figure_7B5<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="SIRD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = 0.007'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="40-64 years" & prev_year_race_age$group=="SIRD",]))

#Collegue or Higher
#TYPE 2 DIABETES

Sup_Figure_7C1<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Overall" & prev_year_race_age$age_cat=="≥65 years",], aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("≥65 years")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 50))+
  geom_text(aes(x = 2,y = 48,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Overall" & prev_year_race_age$age_cat=="≥65 years",]))

#MARD
Sup_Figure_7C2<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="MARD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 15))+
  geom_text(aes(x = 2,y = 13,label = 'p value = 0.217'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="MARD",]))

#MOD
Sup_Figure_7C3<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="MOD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 15))+
  geom_text(aes(x = 2,y = 13,label = 'p value = 0.003'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="MOD",]))
#SIID
Sup_Figure_7C4<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="SIID",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIID")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 20))+
  geom_text(aes(x = 2,y = 18,label = 'p value = 0.0108'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="SIID",]))

#SIRD
Sup_Figure_7C5<-ggplot(prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="SIRD",], 
                       aes(x=YEAR_2, y=prop, group=RACE_REC, colour=RACE_REC)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Ethnicity")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 20))+
  geom_text(aes(x = 2,y = 18,label = 'p value = 0.065'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ RACE_REC * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_race_age[prev_year_race_age$subgroup=="Subgroups" & prev_year_race_age$age_cat=="≥65 years" & prev_year_race_age$group=="SIRD",]))

Sup_Figure_7A2_LOW<-ggarrange(Sup_Figure_7A2,Sup_Figure_7A3,Sup_Figure_7A4,Sup_Figure_7A5,ncol = 2,nrow = 2,common.legend = T)
Sup_Figure_7A<-ggarrange(Sup_Figure_7A1,Sup_Figure_7A2_LOW,ncol = 1,nrow = 2,common.legend = T)

Sup_Figure_7B2_LOW<-ggarrange(Sup_Figure_7B2,Sup_Figure_7B3,Sup_Figure_7B4,Sup_Figure_7B5,ncol = 2,nrow = 2,common.legend = T)
Sup_Figure_7B<-ggarrange(Sup_Figure_7B1,Sup_Figure_7B2_LOW,ncol = 1,nrow = 2,common.legend = T)

Sup_Figure_7C2_LOW<-ggarrange(Sup_Figure_7C2,Sup_Figure_7C3,Sup_Figure_7C4,Sup_Figure_7C5,ncol = 2,nrow = 2,common.legend = T)
Sup_Figure_7C<-ggarrange(Sup_Figure_7C1,Sup_Figure_7C2_LOW,ncol = 1,nrow = 2,common.legend = T)

Sup_Figure_7<-ggarrange(Sup_Figure_7A,Sup_Figure_7B,Sup_Figure_7C,ncol = 3,nrow = 1,common.legend = T)

ggsave(Sup_Figure_7,
       filename = "Supp_Fig_7.jpg", 
       width = 50, 
       height = 28,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)


#----Age categorization (Supplementary Figure 4)-----

#Age Stratification
Figure4A<-ggplot(prev_year_age[prev_year_age$subgroup=="Overall",], 
                 aes(x=YEAR_2, y=prop, group=age_cat, colour=age_cat)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Age \nCategories")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 28,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ as.numeric(age_cat) + Error(YEAR_2), data = prev_year_age[prev_year_age$subgroup=="Overall",]))

Figure4B<-ggplot(prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="MARD",], 
                 aes(x=YEAR_2, y=prop, group=age_cat, colour=age_cat)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Age \nCategories")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 15))+
  geom_text(aes(x = 2,y = 13,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ age_cat * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="MARD",]))

Figure4C<-ggplot(prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="MOD",], 
                 aes(x=YEAR_2, y=prop, group=age_cat, colour=age_cat)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Age \nCategories")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 15))+
  geom_text(aes(x = 2,y = 13,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ age_cat * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="MOD",]))


Figure4D<-ggplot(prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="SIRD",], 
                 aes(x=YEAR_2, y=prop, group=age_cat, colour=age_cat)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Age \nCategories")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ age_cat * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="SIRD",]))


Figure4E<-ggplot(prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="SIID",], 
                 aes(x=YEAR_2, y=prop, group=age_cat, colour=age_cat)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Age \nCategories")+ggtitle("SIDD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ age_cat * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_age[prev_year_age$subgroup=="Subgroups" & prev_year_age$group=="SIID",]))


Figure4_Down<-ggpubr::ggarrange(Figure4B,Figure4C,Figure4D,Figure4E,nrow = 2, ncol = 2,common.legend = T)
Figure4<-ggpubr::ggarrange(Figure4A,Figure4_Down,nrow = 1, ncol = 2,labels = c("A","B"),common.legend = T)

ggsave(Figure4,
       filename = "Supp_Fig_4.jpg", 
       width = 35, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)


#----Obesity stratification (Supplementary Figure 9)----

DM2_prev_obesity<-svyby(~DM2_final,~YEAR_2+BMI_CAT, NHANES_subset, svymean) %>%   mutate(DM2_final = round(DM2_final*100, digits=1), se=round(se*100, digits=1))
DM2_MARD_obesity<-svyby(~ClusterMARD, ~YEAR_2+BMI_CAT, NHANES_subset1, svymean) %>%   mutate(ClusterMARD = round(ClusterMARD*100, digits=1), se=round(se*100, digits=1))
DM2_MOD_obesity<-svyby(~ClusterMOD, ~YEAR_2+BMI_CAT , NHANES_subset1, svymean)%>%   mutate(ClusterMOD = round(ClusterMOD*100, digits=1), se=round(se*100, digits=1))
DM2_SIRD_obesity<-svyby(~ClusterSIRD,~YEAR_2+BMI_CAT, NHANES_subset1, svymean) %>%  mutate(ClusterSIRD = round(ClusterSIRD*100, digits=1), se=round(se*100, digits=1))
DM2_SIID_obesity<-svyby(~ClusterSIID, ~YEAR_2+BMI_CAT , NHANES_subset1, svymean)%>%  mutate(ClusterSIID = round(ClusterSIID*100, digits=1), se=round(se*100, digits=1))

names(DM2_prev_obesity)[names(DM2_prev_obesity) == "DM2_final"] <- "prop"
names(DM2_MARD_obesity)[names(DM2_MARD_obesity) == "ClusterMARD"] <- "prop"
names(DM2_MOD_obesity)[names(DM2_MOD_obesity) == "ClusterMOD"] <- "prop"
names(DM2_SIRD_obesity)[names(DM2_SIRD_obesity) == "ClusterSIRD"] <- "prop"
names(DM2_SIID_obesity)[names(DM2_SIID_obesity) == "ClusterSIID"] <- "prop"

DM2_prev_obesity$group<-"Type 2 Diabetes"
DM2_MARD_obesity$group<-"MARD"
DM2_MOD_obesity$group<-"MOD"
DM2_SIRD_obesity$group<-"SIRD"
DM2_SIID_obesity$group<-"SIID"

prev_year_obesity<-rbind(DM2_prev_obesity,DM2_MARD_obesity,DM2_MOD_obesity,DM2_SIRD_obesity,DM2_SIID_obesity)
prev_year_obesity$subgroup<-ifelse(prev_year_obesity$group=="Type 2 Diabetes", 0, 1)
prev_year_obesity$subgroup<-factor(prev_year_obesity$subgroup, labels = c("Overall","Subgroups"))

prev_year_obesity$low_limit<-prev_year_obesity$prop-prev_year_obesity$se
prev_year_obesity$upp_limit<-prev_year_obesity$prop+prev_year_obesity$se
prev_year_obesity$IC<-paste0(str_pad(prev_year_obesity$prop, 1,pad = "0"),
                             c(" "),
                             c("("),
                             str_pad(prev_year_obesity$low_limit, 1,pad = "0"),
                             c("-"),
                             str_pad(prev_year_obesity$upp_limit,1, pad="0"),
                             c(")"))

write.csv(prev_year_obesity, "prev_year_obesity.csv")

#BMI Stratification

Figure9A<-ggplot(prev_year_obesity[prev_year_obesity$subgroup=="Overall",], aes(x=YEAR_2, y=prop, group=BMI_CAT, colour=BMI_CAT)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="BMI Categories")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 30))+
  geom_text(aes(x = 2,y = 25,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ BMI_CAT * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_obesity[prev_year_obesity$subgroup=="Overall",]))

Figure9B<-ggplot(prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="MARD",], aes(x=YEAR_2, y=prop, group=BMI_CAT, colour=BMI_CAT)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="BMI Categories")+ggtitle("MARD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value = < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ BMI_CAT * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="MARD",]))

Figure9C<-ggplot(prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="MOD",], aes(x=YEAR_2, y=prop, group=BMI_CAT, colour=BMI_CAT)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="BMI Categories")+ggtitle("MOD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0,15))+
  geom_text(aes(x = 2,y = 13,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ BMI_CAT * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="MOD",]))

Figure9D<-ggplot(prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="SIRD",], aes(x=YEAR_2, y=prop, group=BMI_CAT, colour=BMI_CAT)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="BMI Categories")+ggtitle("SIRD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ BMI_CAT * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="SIRD",]))

Figure9E<-ggplot(prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="SIID",], aes(x=YEAR_2, y=prop, group=BMI_CAT, colour=BMI_CAT)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="BMI Categories")+ggtitle("SIDD")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

summary(aov(as.numeric(prop) ~ BMI_CAT * as.numeric(YEAR_2) + Error(YEAR_2), data = prev_year_obesity[prev_year_obesity$subgroup=="Subgroups" &prev_year_obesity$group=="SIID",]))

Figure9_Down<-ggpubr::ggarrange(Figure9B,Figure9C,Figure9D,Figure9E,nrow = 2, ncol = 2,common.legend = T)
Figure9<-ggpubr::ggarrange(Figure9A,Figure9_Down,nrow = 1, ncol = 2,labels = c("A","B"),common.legend = T)

ggsave(Figure9,
       filename = "Supp_Fig_5.jpg", 
       width = 35, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)


#----Recently diagnosed diabetes (Supplementary Figure 10)----

DM2_prev_recent<-svyby(~RECENT_DIABETES,~YEAR_2, NHANES_subset, svymean) %>%   mutate(RECENT_DIABETES1 = round(RECENT_DIABETES1*100, digits=1), 
                                                                                      se.RECENT_DIABETES1=round(se.RECENT_DIABETES1*100, digits=1),
                                                                                      RECENT_DIABETES2 = round(RECENT_DIABETES2*100, digits=1), 
                                                                                      se.RECENT_DIABETES2=round(se.RECENT_DIABETES2*100, digits=1),
                                                                                      RECENT_DIABETES3 = round(RECENT_DIABETES3*100, digits=1), 
                                                                                      se.RECENT_DIABETES3=round(se.RECENT_DIABETES3*100, digits=1))

DM2_MARD_recent<-svyby(~RECENT_MARD, ~YEAR_2, NHANES_subset1, svymean) %>%   mutate(RECENT_MARD1 = round(RECENT_MARD1*100, digits=1), 
                                                                                    se.RECENT_MARD1=round(se.RECENT_MARD1*100, digits=1),
                                                                                    RECENT_MARD2 = round(RECENT_MARD2*100, digits=1), 
                                                                                    se.RECENT_MARD2=round(se.RECENT_MARD2*100, digits=1),
                                                                                    RECENT_MARD3 = round(RECENT_MARD3*100, digits=1), 
                                                                                    se.RECENT_MARD3=round(se.RECENT_MARD3*100, digits=1),
                                                                                    RECENT_MARD4 = round(RECENT_MARD4*100, digits=1), 
                                                                                    se.RECENT_MARD4=round(se.RECENT_MARD4*100, digits=1),
                                                                                    RECENT_MARD5 = round(RECENT_MARD5*100, digits=1), 
                                                                                    se.RECENT_MARD5=round(se.RECENT_MARD5*100, digits=1))

DM2_MOD_recent<-svyby(~RECENT_MOD, ~YEAR_2, NHANES_subset1, svymean) %>%   mutate(RECENT_MOD = round(RECENT_MOD1*100, digits=1), 
                                                                                  se.RECENT_MOD1=round(se.RECENT_MOD1*100, digits=1),
                                                                                  RECENT_MOD2 = round(RECENT_MOD2*100, digits=1), 
                                                                                  se.RECENT_MOD2=round(se.RECENT_MOD2*100, digits=1),
                                                                                  RECENT_MOD3 = round(RECENT_MOD3*100, digits=1), 
                                                                                  se.RECENT_MOD3=round(se.RECENT_MOD3*100, digits=1),
                                                                                  RECENT_MOD4 = round(RECENT_MOD4*100, digits=1), 
                                                                                  se.RECENT_MOD4=round(se.RECENT_MOD4*100, digits=1),
                                                                                  RECENT_MOD5 = round(RECENT_MOD5*100, digits=1), 
                                                                                  se.RECENT_MOD5=round(se.RECENT_MOD5*100, digits=1))

DM2_SIRD_recent<-svyby(~RECENT_SIRD, ~YEAR_2, NHANES_subset1, svymean) %>%   mutate(RECENT_SIRD = round(RECENT_SIRD1*100, digits=1), 
                                                                                    se.RECENT_SIRD1=round(se.RECENT_SIRD1*100, digits=1),
                                                                                    RECENT_SIRD2 = round(RECENT_SIRD2*100, digits=1), 
                                                                                    se.RECENT_SIRD2=round(se.RECENT_SIRD2*100, digits=1),
                                                                                    RECENT_SIRD3 = round(RECENT_SIRD3*100, digits=1), 
                                                                                    se.RECENT_SIRD3=round(se.RECENT_SIRD3*100, digits=1),
                                                                                    RECENT_SIRD4 = round(RECENT_SIRD4*100, digits=1), 
                                                                                    se.RECENT_SIRD4=round(se.RECENT_SIRD4*100, digits=1),
                                                                                    RECENT_SIRD5 = round(RECENT_SIRD5*100, digits=1), 
                                                                                    se.RECENT_SIRD5=round(se.RECENT_SIRD5*100, digits=1))

DM2_SIID_recent<-svyby(~RECENT_SIID, ~YEAR_2, NHANES_subset1, svymean) %>%   mutate(RECENT_SIID = round(RECENT_SIID1*100, digits=1), 
                                                                                    se.RECENT_SIID1=round(se.RECENT_SIID1*100, digits=1),
                                                                                    RECENT_SIID2 = round(RECENT_SIID2*100, digits=1), 
                                                                                    se.RECENT_SIID2=round(se.RECENT_SIID2*100, digits=1),
                                                                                    RECENT_SIID3 = round(RECENT_SIID3*100, digits=1), 
                                                                                    se.RECENT_SIID3=round(se.RECENT_SIID3*100, digits=1),
                                                                                    RECENT_SIID4 = round(RECENT_SIID4*100, digits=1), 
                                                                                    se.RECENT_SIID4=round(se.RECENT_SIID4*100, digits=1),
                                                                                    RECENT_SIID5 = round(RECENT_SIID5*100, digits=1), 
                                                                                    se.RECENT_SIID5=round(se.RECENT_SIID5*100, digits=1))

Figure10A<-ggplot2::ggplot(DM2_prev_recent) + 
  geom_line(aes(x=YEAR_2, y=RECENT_DIABETES2,group=1, col="Recently Diagnosed Diabetes"),size=2)+
  geom_line(aes(x=YEAR_2, y=RECENT_DIABETES3,group=1, col="≥5 Years from Diabetes Diagnosis"),size=2)+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_DIABETES2+se.RECENT_DIABETES2, ymax=RECENT_DIABETES2-se.RECENT_DIABETES2), width=.2,
                position=position_dodge(0.01))+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_DIABETES3+se.RECENT_DIABETES3, ymax=RECENT_DIABETES3-se.RECENT_DIABETES3), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  labs(colour="")+ggtitle("Overall")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_color_jama()+ 
  scale_y_continuous(limits = c(0, 10))+
  geom_text(aes(x = 2,y = 8,label = 'p value <0.001'),show.legend = F,check_overlap = TRUE)

Figure10B<-ggplot2::ggplot(DM2_MARD_recent) + 
  geom_line(aes(x=YEAR_2, y=RECENT_MARD4,group=1, col="Recently Diagnosed"),size=2)+
  geom_line(aes(x=YEAR_2, y=RECENT_MARD5,group=1, col="≥5 Years from Diabetes Diagnosis"),size=2)+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_MARD4+se.RECENT_MARD4, ymax=RECENT_MARD4-se.RECENT_MARD4), width=.2,
                position=position_dodge(0.01))+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_MARD5+se.RECENT_MARD5, ymax=RECENT_MARD5-se.RECENT_MARD5), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+ggtitle("MARD")+
  labs(colour="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_color_jama()+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

Figure10C<-ggplot2::ggplot(DM2_MOD_recent) + 
  geom_line(aes(x=YEAR_2, y=RECENT_MOD4,group=1, col="Recently Diagnosed"),size=2)+
  geom_line(aes(x=YEAR_2, y=RECENT_MOD5,group=1, col="≥5 Years from Diabetes Diagnosis"),size=2)+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_MOD4+se.RECENT_MOD4, ymax=RECENT_MOD4-se.RECENT_MOD4), width=.2,
                position=position_dodge(0.01))+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_MOD5+se.RECENT_MOD5, ymax=RECENT_MOD5-se.RECENT_MOD5), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+ggtitle("MOD")+
  labs(colour="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_color_jama()+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

Figure10D<-ggplot2::ggplot(DM2_SIRD_recent) + 
  geom_line(aes(x=YEAR_2, y=RECENT_SIRD4,group=1, col="Recently Diagnosed"),size=2)+
  geom_line(aes(x=YEAR_2, y=RECENT_SIRD5,group=1, col="≥5 Years from Diabetes Diagnosis"),size=2)+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_SIRD4+se.RECENT_SIRD4, ymax=RECENT_SIRD4-se.RECENT_SIRD4), width=.2,
                position=position_dodge(0.01))+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_SIRD5+se.RECENT_SIRD5, ymax=RECENT_SIRD5-se.RECENT_SIRD5), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+ggtitle("SIRD")+
  labs(colour="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_color_jama()+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)

Figure10E<-ggplot2::ggplot(DM2_SIID_recent) + 
  geom_line(aes(x=YEAR_2, y=RECENT_SIID4,group=1, col="Recently Diagnosed"),size=2)+
  geom_line(aes(x=YEAR_2, y=RECENT_SIID5,group=1, col="≥5 Years from Diabetes Diagnosis"),size=2)+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_SIID4+se.RECENT_SIID4, ymax=RECENT_SIID4-se.RECENT_SIID4), width=.2,
                position=position_dodge(0.01))+
  geom_errorbar(aes(x=YEAR_2,ymin=RECENT_SIID5+se.RECENT_SIID5, ymax=RECENT_SIID5-se.RECENT_SIID5), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+ggtitle("SIDD")+
  labs(colour="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))+
  scale_color_jama()+ 
  scale_y_continuous(limits = c(0, 5))+
  geom_text(aes(x = 2,y = 4.5,label = 'p value < 0.001'),show.legend = F,check_overlap = TRUE)


Figure10_Down<-ggpubr::ggarrange(Figure10B,Figure10C,Figure10D,Figure10E,nrow = 2, ncol = 2,common.legend = T)
Figure10<-ggpubr::ggarrange(Figure10A,Figure10_Down,nrow = 1, ncol = 2,labels = c("A","B"),common.legend = T)

ggsave(Figure10,
       filename = "Supp_Fig_6.jpg", 
       width = 35, 
       height = 15,
       units=c("cm"),
       dpi = 300,
       limitsize = FALSE)


#--- Piecewise reg----

#Logistic models
set.seed(123)
HANES_DM2$YEAR_NUM<-NULL;
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="1999-2000"]<-1
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2001-2002"]<-2
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2003-2004"]<-3
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2005-2006"]<-4
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2007-2008"]<-5
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2009-2010"]<-6
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2011-2012"]<-7
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2013-2014"]<-8
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2015-2016"]<-9
NHANES_DM2$YEAR_NUM[NHANES_DM2$YEAR=="2017-2018"]<-10
NHANES_DM2$YEAR_NUM<-as.numeric(NHANES_DM2$YEAR_NUM)

log1<-glm(DM2_final~YEAR_NUM,  data=NHANES_DM2, family="poisson")
log2<-glm(ClusterMARD~YEAR_NUM, data=NHANES_DM2, family="poisson")
log3<-glm(ClusterMOD~YEAR_NUM, data=NHANES_DM2, family="poisson")
log4<-glm(ClusterSIID~YEAR_NUM, data=NHANES_DM2, family="poisson")
log5<-glm(ClusterSIRD~YEAR_NUM, data=NHANES_DM2, family="poisson")

#Breakpoints identification
set.seed(123)

davies.test(log1, seg.Z=~YEAR_NUM)
fit_seg1<-segmented(log1, seg.Z=~YEAR_NUM, psi=list(YEAR_NUM=c(2)), npsi=5)
print(fit_seg1, include.psi=T)
davies.test(fit_seg1, seg.Z=~YEAR_NUM)
plot(fit_seg1, conf.level=.95, is=TRUE, isV=TRUE, col=1, shade = TRUE, col.shade=1)
slope(fit_seg1)

davies.test(log2, seg.Z=~YEAR_NUM) #Log2 (MARD) no change in slope
fit_seg2<-segmented(log2, seg.Z=~YEAR_NUM, psi=list(YEAR_NUM=c(2)), npsi=5)
print(fit_seg2, include.psi=T)
davies.test(fit_seg2, seg.Z=~YEAR_NUM)
plot(fit_seg2, conf.level=.9, is=TRUE, isV=TRUE, col=2, shade = TRUE, col.shade=2)

davies.test(log3, seg.Z=~YEAR_NUM) #MOD change in slope
fit_seg3<-segmented(log3, seg.Z=~YEAR_NUM, psi=list(YEAR_NUM=c(2,3)), npsi=5)
print(fit_seg3, include.psi=T)
davies.test(fit_seg3, seg.Z=~YEAR_NUM)
plot(fit_seg3, conf.level=.9, is=TRUE, isV=TRUE, col=3, shade = TRUE, col.shade=3)

davies.test(log4, seg.Z=~YEAR_2)#log (SIID) No change in slope
fit_seg4<-segmented(log4, seg.Z=~YEAR_2, psi=psis, npsi=5)
print(fit_seg4, include.psi=T)
davies.test(fit_seg4, seg.Z=~YEAR_2)
plot(fit_seg4, conf.level=.9, is=TRUE, isV=TRUE, col=4, shade = TRUE, col.shade=4)

davies.test(log5, seg.Z=~YEAR_2) #SIRD change in slope
fit_seg5<-segmented(log5, seg.Z=~YEAR_2, psi=list(YEAR_2=c(2,3,4)), npsi=5)
print(fit_seg5, include.psi=T)
davies.test(fit_seg5, seg.Z=~YEAR_2)
plot(fit_seg5, conf.level=.9, is=TRUE, isV=TRUE, col=5, shade = TRUE, col.shade=5)



#PRevious dataset

NHANES_DM2<-readr::read_rds("NHANES_DM2_FINAL.rds")


