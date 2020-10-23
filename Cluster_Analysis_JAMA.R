#Trends in the prevalence of diabetes subgroups in U.S adults: A data-driven cluster analysis in NHANES from 1998 to 2018 
#Analysis: Neftali Eduardo Antonio-Villa (nefoantonio@hotmail.com); Luisa Fernandez-Chirino (); Omar Yaxmehen Bello-Chavolla (oyaxbell@yahoo.com.mx)
#
#Diclosure: 

#----Library loading----

library(haven); 
library(tidyverse); library(dplyr); 
library(survival)
library(ggthemes)
library(ggsci)
library(survey)
library(segmented)

#----Dataset managment-----
setwd("/Users/nefoantonio/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/PROYECTOS/T2D CLUSTERS")
NHANES_DM2<-readr::read_rds("NHANES_DM2.rds")

table(NHANES_DM2$DM2_final)
nrow(NHANES_DM2)

NHANES_DM2$DM2_final<-na.tools::na.replace(NHANES_DM2$DM2_final,0)
NHANES_DM2$AGE_18<-NULL;NHANES_DM2$AGE_18[NHANES_DM2$AGE>=18]<-1;NHANES_DM2$AGE_18[NHANES_DM2$AGE<18]<-0
NHANES_DM2$BMI_2<-NHANES_DM2$Weight/((NHANES_DM2$Height/100)^2)

NHANES_DM2<-NHANES_DM2%>%mutate(one = 1)

d1<-dummies::dummy(NHANES_DM2$Cluster)
NHANES_DM2<-cbind(NHANES_DM2,d1)

table(NHANES_DM2$DM2_final, NHANES_DM2$AGE_18)
table(NHANES_DM2$Cluster, NHANES_DM2$AGE_18)


#Recoding of variables
NHANES_DM2$YEAR_2<-NULL;
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="1999-2000"]<-1
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2001-2002"]<-1
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2003-2004"]<-2
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2005-2006"]<-2
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2007-2008"]<-3
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2009-2010"]<-3
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2011-2012"]<-4
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2013-2014"]<-4
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2015-2016"]<-5
NHANES_DM2$YEAR_2[NHANES_DM2$YEAR=="2017-2018"]<-5
NHANES_DM2$YEAR_2<-factor(NHANES_DM2$YEAR_2, labels = c("1999-2002","2003-2006","2007-2010","2011-2014","2015-2018"))


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

#----Figure 1-----
Figure1A<-ggplot(prev_year, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+facet_wrap(~subgroup, scales = "free")

Figure1B<-ggplot(prev_year_race, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~subgroup+RACE_REC, nrow=2, scales = "free")

Figure1<-ggpubr::ggarrange(Figure1A,Figure1B,nrow = 1, ncol = 2,labels = c("A","B"),common.legend = T)

ggsave(Figure1,
       filename = "Figure1.jpg", 
       width = 65, 
       height = 15,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

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

#----Obesity stratification----


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

#----Figure 2----

Figure2A<-ggplot(prev_year_age, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~age_cat+subgroup, nrow=3, scales = "free")


Figure2B<-ggplot(prev_year_obesity, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~BMI_CAT+subgroup, nrow=3, scales = "free")

Figure2<-ggpubr::ggarrange(Figure2A,Figure2B,nrow = 1, ncol = 2,labels = c("A","B"),common.legend = T)

ggsave(Figure2,
       filename = "Figure2.jpg", 
       width = 65, 
       height = 20,
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

ggplot(prev_year_race_age_obesity, aes(x=YEAR_2, y=prop, group=group, colour=group)) + 
  geom_line(size=1.5) +
  geom_point(size=2)+
  geom_errorbar(aes(ymin=prop+se, ymax=prop-se), width=.2,
                position=position_dodge(0.01))+
  theme_hc()+
  ylab("Prevalence, (%)")+
  xlab("NHANES cycle")+
  scale_color_jama()+
  labs(colour="Groups")+
  facet_wrap(~RACE_REC+age_cat+subgroup+BMI_CAT, nrow=3, scales = "free")


#----Trend analysis-----
library(quantreg)

summay(rq(prev~year, data=data))
#Overall diabetes and subgroups trends
summary(glm(prop~as.numeric(YEAR_2), data=DM2_prev))
summary(glm(prop~as.numeric(YEAR_2), data=DM2_MARD))
summary(glm(prop~as.numeric(YEAR_2), data=DM2_MOD))
summary(glm(prop~as.numeric(YEAR_2), data=DM2_SIRD))
summary(glm(prop~as.numeric(YEAR_2), data=DM2_SIID))

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
summary(rq(prop~as.numeric(YEAR_2), data=prev_year_race, subset = (subgroup=="Subgroups" & RACE_REC=="Mexican American" & group=="MOD")))
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

#----Obesity weighted prevalence-----

NHANES_DM2$OBESITY_2<-NULL;
NHANES_DM2$OBESITY_2[NHANES_DM2$BMI_2<30]<-0
NHANES_DM2$OBESITY_2[NHANES_DM2$BMI_2>=30]<-1

NHANES_DM2_ob<-NHANES_DM2%>%filter(BMI_2>=0)

NHANES_all_ob <- svydesign(data=NHANES_DM2_ob, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
NHANES_subset_ob <- subset(NHANES_all_ob, BMI_2>=0)
Ob_prev<-svyby(~OBESITY_2, ~YEAR_2, NHANES_subset_ob, svymean) %>%   mutate(OBESITY_2 = round(OBESITY_2*100, digits=1), se=round(se*100, digits=1))
View(Ob_prev)


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
