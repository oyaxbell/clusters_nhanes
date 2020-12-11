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


#----Dataset managment-----
setwd("/Users/nefoantonio/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/PROYECTOS/T2D CLUSTERS")
NHANES_DM2<-readr::read_rds("NHANES_DM2_III_2019_2.rds")

#Descriptive Analysis
nrow(NHANES_DM2)
table(NHANES_DM2$AGE>=18)
NHANES_DM2$Glucose[NHANES_DM2$Glucose == 88888] <- NA
NHANES_DM2$Glycohemoglobin[NHANES_DM2$Glycohemoglobin == 8888] <- NA
NHANES_DM2$DM2_final<-as.numeric((NHANES_DM2$Glycohemoglobin>=6.5 |NHANES_DM2$Glucose>=126 | NHANES_DM2$Dr_diabetes==1))
NHANES_DM2$DM2_final<-na.tools::na.replace(NHANES_DM2$DM2_final,0)

table(NHANES_DM2$DM2_final,NHANES_DM2$AGE>=18)
sum(!is.na(NHANES_DM2$Cluster))
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

nhanes3<-read_csv("/Users/nefoantonio/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/PROYECTOS/T2D CLUSTERS/nhanes3.csv")
nhanes3_edu<-nhanes3%>%dplyr::select(SEQN,SDPPHASE.x,HFA8R)
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

edu_n<-read_csv("Letter CLUSTERS/Datasets CSV/educ.csv")
edu_n$seqn<-paste0(str_pad(edu_n$SEQN, 1,pad = "0"),c("-4"))
edu_n$EDU_REC<-edu_n$Schooling

edu_n_2<-edu_n%>%dplyr::select(seqn,EDU_REC)
edu_final<-rbind(nhanes3_edu_2,edu_n_2)

NHANES_DM2<-NHANES_DM2%>%left_join(edu_final,by="seqn")
NHANES_DM2$EDU_REC_2<-NULL
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==1]<-1
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==2]<-2
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==3]<-2
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==4]<-3
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==5]<-3
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==7]<-1
NHANES_DM2$EDU_REC_2[NHANES_DM2$EDU_REC==9]<-1

NHANES_DM2$EDU_REC_2<-factor(NHANES_DM2$EDU_REC_2, labels = c("Primary School or No-Education", "Secondary, High-School or AA Degre", "Collegue or Higher"))


#----Diabetes weighted subset----

NHANES_all <- svydesign(data=NHANES_DM2, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
NHANES_subset <- subset(NHANES_all, AGE_18==1)
NHANES_subset1<- subset(NHANES_all, AGE_18==1 & !is.na(NHANES_DM2$Insulin))

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
DM2_SIID$group<-"SIID"

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
DM2_SIID_sex$group<-"SIID"

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

write.csv(prev_year_sex, "prev_year_sex.csv")
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

write.csv(prev_year_race, "prev_year_race.csv")

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

write.csv(prev_year_age, "prev_year_age.csv")

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

write.csv(prev_year_edu, "prev_year_edu.csv")
#----Figure 1-----
Figure1A<-ggplot(prev_year, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  scale_color_jama()+
  labs(colour="Groups")+facet_wrap(~subgroup, scales = "free")

Figure1B<-ggplot(prev_year_sex, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~subgroup+GENDER_REC, nrow=2, scales = "free")

Figure2C<-ggplot(prev_year_race, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~subgroup+RACE_REC, nrow=2, scales = "free")

Figure2D<-ggplot(prev_year_age, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~subgroup+age_cat, nrow=2, scales = "free")

Figure2E<-ggplot(prev_year_edu, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~subgroup+EDU_REC_2, nrow=2, scales = "free")


FIGURE_1_VERTICAL<-ggpubr::ggarrange(Figure1A,Figure1B,nrow = 2, ncol = 1,labels = c("A","B"),common.legend = T)
FIGURE_2_VERTICAL<-ggpubr::ggarrange(Figure2C,Figure2E,Figure2D,nrow = 3, ncol = 1,labels = c("A","B","C"),common.legend = T)

ggsave(FIGURE_1_VERTICAL,
       filename = "Figure1.jpg", 
       width = 40, 
       height = 25,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

ggsave(FIGURE_2_VERTICAL,
       filename = "Figure2.jpg", 
       width = 40, 
       height = 40,
       units=c("cm"),
       dpi = 500,
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



#----Education and etnicity stratification (Sup Fig 3)----
 
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
DM2_SIID_race_edu$group<-"SIID"

prev_year_race_edu<-rbind(DM2_prev_race_edu,DM2_MARD_race_edu,DM2_MOD_race_edu,DM2_SIRD_race_edu,DM2_SIID_race_edu)
prev_year_race_edu$subgroup<-ifelse(prev_year_race_edu$group=="Type 2 Diabetes", 0, 1)
prev_year_race_edu$subgroup<-factor(prev_year_race_edu$subgroup, labels = c("Overall","Subgroups"))

Supp_Fig_X<-ggplot(prev_year_race_edu, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~RACE_REC+EDU_REC_2+subgroup, nrow=3, scales = "free")

ggsave(Supp_Fig_X,
       filename = "Supp_Fig_10.jpg", 
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

#----Obesity stratification (Sup Fig 4)----

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

Sup_Figure8<-ggplot(prev_year_obesity, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~subgroup+BMI_CAT, nrow=2, scales = "free")

ggsave(Sup_Figure8,
       filename = "Supp_Fig_8.jpg", 
       width = 60, 
       height = 15,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)



#----Ethnicity and age stratification-----

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
DM2_SIID_race_age$group<-"SIID"

prev_year_race_age<-rbind(DM2_prev_race_age,DM2_MARD_race_age,DM2_MOD_race_age,DM2_SIRD_race_age,DM2_SIID_race_age)
prev_year_race_age$subgroup<-ifelse(prev_year_race_age$group=="Type 2 Diabetes", 0, 1)
prev_year_race_age$subgroup<-factor(prev_year_race_age$subgroup, labels = c("Overall","Subgroups"))


Supp_Fig_7<-ggplot(prev_year_race_age, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~RACE_REC+age_cat+subgroup, nrow=3, scales = "free")

ggsave(Supp_Fig_7,
       filename = "Supp_Fig_7.jpg", 
       width = 75, 
       height = 20,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

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


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Primary School or No-Education")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Primary School or No-Education")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Secondary, High-School or AA Degre")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Secondary, High-School or AA Degre")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Collegue or Higher")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Overall" & EDU_REC_2=="Collegue or Higher")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Primary School or No-Education" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="MARD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="MARD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="MOD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="MOD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_edu, subset = (subgroup=="Subgroups" & EDU_REC_2=="Collegue or Higher" & group=="SIID")))


#Education and etinicity stratification
table(prev_year_race_edu$RACE_REC)
summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Primary School or No-Education" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Primary School or No-Education" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIRD")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Collegue or Higher" & group=="SIRD")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & EDU_REC_2=="Collegue or Higher" & group=="SIRD")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Primary School or No-Education" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Primary School or No-Education" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Collegue or Higher" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic White" & EDU_REC_2=="Collegue or Higher" & group=="SIID")))


summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Primary School or No-Education" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Primary School or No-Education" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Secondary, High-School or AA Degre" & group=="SIID")))

summary(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Collegue or Higher" & group=="SIID")))
confint(glm(prop~as.numeric(YEAR_2), data=prev_year_race_edu, subset = (subgroup=="Subgroups" & RACE_REC=="Non-Hispanic Black" & EDU_REC_2=="Collegue or Higher" & group=="SIID")))





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

#----Proportion of diabetes in te studied period (Sup Fig 5)-----

NHANES_DM2_Prev<-NHANES_DM2%>%dplyr::select(Cluster,YEAR)%>%na.omit()

NHANES_DM2_Prev_FIG<-NHANES_DM2_Prev%>%
  group_by(YEAR,Cluster)%>%
  summarise(freq_clus=n())%>%
  mutate(prop_clus_year=(freq_clus/sum(freq_clus))*100)%>%
  na.omit()

Supp_Fig_5<-ggplot(NHANES_DM2_Prev_FIG, aes(y=prop_clus_year, x=Cluster,fill=Cluster)) +
  geom_bar(stat="identity", width=0.8, position=position_dodge(width=0.8))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme_hc()+
  ylab("Absolute Proportion (%)")+
  geom_text(position = position_dodge(width= 0.8),aes(label=round(prop_clus_year*1,2)), vjust=1.6,color="white", size=3)+
  xlab("NHANES Cycle")+
  scale_fill_jama()+
  labs(fill="Diabetes Clusters")+
  facet_wrap(~YEAR, nrow = 3)

ggsave(Supp_Fig_5,
       filename = "Supp_Fig_5.jpg", 
       width = 45, 
       height = 20,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

#----Differences of diabetes (Sup Fig 6)-----

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

Sup_Fig_3A<-ggplot(data, aes(x="", y=`Absolute Frecuency`, fill=`Diabetes Cluster`)) +
  geom_bar(stat="identity", width=2) +
  coord_polar("y", start=0)+
  theme_void()+
  ylab("")+
  theme(legend.position="bottom") +
  geom_text(aes(y = ypos, label = abs_freq_2), color = "white", size=6) +
  scale_fill_jama()

#Hba1c
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
  ylab("HOMA2-2B")+
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
       filename = "Supp_Fig_6.jpg", 
       width = 80, 
       height = 30,
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

data2<-as.data.frame(prop.table(table(NHANES_DM2_Cluster[NHANES_DM2_Cluster$RACE_REC=="Non-Hispanic White",]$Cluster))*100)
names(data2)<-c("Diabetes Cluster", "Absolute Frecuency")

data2 <- data2 %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
data2$abs_freq_2<-as.numeric(format(round(data2$`Absolute Frecuency`, 2), nsmall = 2))

data3<-as.data.frame(prop.table(table(NHANES_DM2_Cluster[NHANES_DM2_Cluster$RACE_REC=="Non-Hispanic Black",]$Cluster))*100)
names(data3)<-c("Diabetes Cluster", "Absolute Frecuency")

data3 <- data3 %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
data3$abs_freq_2<-as.numeric(format(round(data3$`Absolute Frecuency`, 2), nsmall = 2))

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

Sup_Fig_4_L<-ggarrange(Sup_Fig_4A.1,Sup_Fig_4A.2,Sup_Fig_4A.3,ncol = 3, nrow = 1, labels = c("A"))

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
  ylab("HOMA2-2B")+
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
Sup_Fig_7<-ggarrange(Sup_Fig_4_L,Sup_Fig_4_R, nrow = 1,ncol = 2, labels = c("A"), common.legend = T)

ggsave(Sup_Fig_7,
       filename = "Supp_Fig_7.jpg", 
       width = 80, 
       height = 35,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

#----SAID stratification in NHANES III (Sup Fig 7)-----

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
DM2_SIID$group<-"SIID"
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
  labs(colour="Groups")+facet_wrap(~subgroup, scales = "free")


NHANES_DM2_Prev_III<-nhanes_3%>%dplyr::select(Cluster_REC,YEAR.x)%>%na.omit()

NHANES_DM2_Prev_FIG<-NHANES_DM2_Prev_III%>%
  group_by(YEAR.x,Cluster_REC)%>%
  summarise(freq_clus=n())%>%
  mutate(prop_clus_year=(freq_clus/sum(freq_clus))*100)%>%
  na.omit()

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
#Frecuncy of clusters
data<-as.data.frame(prop.table(table(NHANES_DM2_Cluster_III$Cluster_REC))*100)
names(data)<-c("Diabetes Cluster", "Absolute Frecuency")

data <- data %>% 
  arrange(desc(`Diabetes Cluster`)) %>%
  mutate(prop = `Absolute Frecuency` / sum(data$`Absolute Frecuency`) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

data$abs_freq_2<-as.numeric(format(round(data$`Absolute Frecuency`, 2), nsmall = 2))

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
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#HOMA2-IR
Sup_Fig_6E<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=HOMA2_IR, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-IR")+
  scale_fill_jama()+
  scale_y_log10()+  
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#HOMA2-B
Sup_Fig_6F<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=HOMA2_B, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("HOMA2-2B")+
  scale_fill_jama()+
  scale_y_log10()+  
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#Age at diagnosis
Sup_Fig_6G<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=Age_diabetes, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Age at Diabetes Diagnosis
(years)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")

#Body Mass Index
Sup_Fig_6H<-ggplot(NHANES_DM2_Cluster_III, aes(x=Cluster_REC,y=BMI, fill=Cluster_REC))+
  geom_boxplot()+
  theme_hc()+
  xlab("Diabetes Cluster")+
  ylab("Body Mass Index (kg/m2)")+
  scale_fill_jama()+
  stat_compare_means(method= "wilcox.test",label = "p.signif",comparisons = my_comparisons)+
  theme(legend.position="none")


Sup_Fig_6_Upp_R<-ggarrange(Supp_Fig_6D,Sup_Fig_6E,Sup_Fig_6F,ncol = 3, nrow = 1, labels = c("D","E","F"))
Sup_Fig_6_Down_R<-ggarrange(Sup_Fig_6G,Sup_Fig_6H,ncol = 2, nrow = 1, labels = c("G","H"))
Sup_Fig_6_LOW_RIGHT<-ggarrange(Sup_Fig_6_Upp_R,Sup_Fig_6_Down_R, nrow = 2,ncol = 1)
Sup_Fig_6_DOWN<-ggarrange(Supp_Fig_6C,Sup_Fig_6_LOW_RIGHT,  labels= c("C"),nrow = 1,ncol = 2)
Sup_Fig_6_UP<-ggarrange(Supp_Fig_6A,Supp_Fig_6B, labels= c("A","B"),nrow = 1,ncol = 2)
Sup_Fig_8<-ggarrange(Sup_Fig_6_UP,Sup_Fig_6_DOWN, nrow = 2,ncol = 1)

ggsave(Sup_Fig_8,
       filename = "Supp_Fig_8.jpg", 
       width = 50, 
       height = 30,
       units=c("cm"),
       dpi = 500,
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


