##################################################################################
##
## Name: seqVA.R
## Description: This file provides an initial analysis for the dataset available 
## for sequence analysis.
## Date: 11/11/19
##
#################################################################################
## HOUSEKEEPING

## set working directory
setwd("//Users/michellemellers/Documents/1_PostDoc/2_VA_Sequence/Test")

## load libraries
## devtools::install_github("ilangurudev/approxmapR",force=TRUE)

library(approxmapR)
library(tidyverse)
library(kableExtra)
library(knitr)
##options(knitr.table.format="html")
library(shiny)
library(lubridate)
##library(plyr)
##library(dplyr)

##################################################################################
## read the data

# evts <- read.csv("./data/eventns_nonphi.csv",stringsAsFactors = F)
# people <- read.csv("./data/companion_nonphi.csv",stringsAsFactors = F)
# ##event2 <- read.csv("./data/eventlg_nonphi.csv",stringsAsFactors = F)
# master <- read.csv("./data/master.csv",stringsAsFactors = F)
# 
# ##################################################################################
# ## Merge the datasets to create a master dataset.
# master <- dplyr::left_join(evts,people,by="randomID")
# ## flag for subset, ie eventsns_nonphi.csv
# ##master <- master %>% mutate(subs = ifelse(randomID %in% evts$randomID,1,0))
# 
# ## create variable for cumulative summary(dt_sum) and if over 90 days(day90)
# master <- master %>% arrange(randomID,order) %>%
#           group_by(randomID) %>% 
#           mutate(dt_sum = cumsum(dt_ID),
#                  days90 = ifelse(dt_sum <= 90,1,0))
# 
# 
# ##write.csv(master,"./data/master.csv",row.names = FALSE)
# 
# ##################################################################################
# ## Answer the scientific question: 
# ## What happens to patients who miss an initial mental illness evaluation?
# ## Dataset, who do not have contact with mental health or in-patient services in the first week.
# ## first subset the data on those who are in 7 days and then see who was not evaluated through 
# ## a missed appointment.
# 
# ## There is a total of 12,067 patients in evts dataset.
# ##master <- master %>% rename("order2"="order")
# 
# ## find who had events within on first week
# evts_sub <- master %>% arrange(randomID,order) %>%
#   group_by(randomID) %>%
#   mutate(servc = ifelse((dt_sum <= 7)&
#      (substr(event,1,2) %in% list("AT","AI"))&
#        (substr(event,4,5) %in% list("MI","MG","MO")),1,0)) %>%
#       filter(servc == 1) %>% distinct(randomID)
# 
# ## every person has an appointment scheduled within the first week
# 
# ## people with mental health appointment scheduled
# evts_sub_MH <- master %>%
#   arrange(randomID,order) %>%
#   group_by(randomID) %>%
#   mutate(servc_MH = ifelse((dt_sum <= 7)&
#                            substr(event,4,5) %in% list("MI","MG","MO"),1,0)) %>%
#   filter(servc_MH == 1) %>% distinct(randomID)
# 
# 
# ## add marker for services in first week not attend
#  master$serv_AT <- ifelse(master$randomID %in% evts_sub$randomID,1,0)
# # ## add marker for scheduled mental illness appointment in first week
#  master$serv_MH <- ifelse(master$randomID %in% evts_sub_MH$randomID,1,0)
# # ## add marker for missing mental health APPT
# master <- master %>%
#          mutate(serv_MISS = ifelse((serv_AT == 0)&(serv_MH == 1),1,0),
#                 serv_MISS3 = case_when(
#                   (serv_AT == 1)&(serv_MH == 1) ~ 1,
#                   (serv_AT == 0)&(serv_MH == 1) ~ 2,
#                   (serv_MH == 0) ~ 0,
#                   TRUE ~ 99
#                 ))
# 
# ## add  variable for discharge quarter and year
# master <- master %>%
#     mutate(disTIME = case_when(
#         (disYR == 2015)&(disQTR == 1) ~ 1,
#         (disYR == 2015)&(disQTR == 2) ~ 2,
#         (disYR == 2015)&(disQTR == 3) ~ 3,
#         (disYR == 2015)&(disQTR == 4) ~ 4,
#         (disYR == 2016)&(disQTR == 1) ~ 5,
#         (disYR == 2016)&(disQTR == 2) ~ 6,
#         (disYR == 2016)&(disQTR == 3) ~ 7,
#         (disYR == 2016)&(disQTR == 4) ~ 8,
#         TRUE ~ 99
#     ),
#       raceEth = case_when(
#         (Race %in% "WHITE")&!(Ethnicity %in% "HISPANIC OR LATINO") ~ "NOT HISPANIC WHITE",
#         (Race %in% "BLACK OR AFRICAN AMERICAN")&!(Ethnicity %in% "HISPANIC OR LATINO") ~ "NOT HISPANIC AFRICAN AMERICAN",
#         (Ethnicity %in% "HISPANIC OR LATINO") ~ "HISPANIC OR LATINO",
#         (Race %in% list("ASIAN","NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER","AMERICAN INDIAN OR ALASKA NATIVE"))&!(Ethnicity %in% "HISPANIC OR LATINO") ~ "NON HISPANIC OTHER",
#         (Race %in% list("","DECLINED TO ANSWER","UNKNOWN BY PATIENT"))&!(Ethnicity %in% "HISPANIC OR LATINO") ~ "UNKNOWN",
#         TRUE ~ "OTHER UNKNOWN"
#       ))
# 
# 
# ## race/ethnicity
# 
# 
# write.csv(master,"./data/master_add.csv",row.names=FALSE)
master <- read.csv("./data/master_add.csv",stringsAsFactors = FALSE)

#############################################################################
## basic descriptive statistics
## comparison to other groups, and with entire dataset

## cross tabs for race
# race <- master %>%
#      filter(order == 1) %>%
#      pivot_longer(c(Race),"Variable","value") %>%
#      group_by(Variable,value,Ethnicity,raceEth) %>%
#      summarize(n=n())

#############################################################################
## cluster the sequences
sel_AIXM <- master %>% filter(days90 == 1) %>%
          mutate(serv_AIXM = ifelse(event == "AI^XM",1,0)) %>%
           filter(serv_AIXM == 1) %>% distinct(randomID)

master$serv_AIXM <- ifelse(master$randomID %in% sel_AIXM$randomID,1,0)  
masterEXC <- master %>% filter((days90 == 1)&(serv_AIXM == 0))


## group if re-entry for mental health
sel_MHR <- masterEXC %>%
            mutate(sel_MHR = ifelse(event %in% c('AI^MI','AI^MG','AI^MO'),1,0)) %>%
          filter(sel_MHR == 1) %>% select(dt_sum,randomID) %>%
          group_by(randomID) %>%
          slice(which.min(dt_sum))

          
masterEXC$sel_MHR <- ifelse((masterEXC$randomID %in% sel_MHR$randomID),1,0)
masterEXC <- dplyr:: left_join(masterEXC,sel_MHR,by=c("randomID")) 
masterEXC$event2 <- ifelse((masterEXC$sel_MHR == 1)&(masterEXC$dt_sum.x <= masterEXC$dt_sum.y),masterEXC$event,NA)
masterEXC$dt_Sum2 <- ifelse((masterEXC$sel_MHR == 1)&(masterEXC$dt_sum.x <= masterEXC$dt_sum.y),masterEXC$dt_sum.x,NA)

## recode 0 to 59 days, 60-90 days, 91+
sel_MHR <- sel_MHR %>%
       mutate(mniQuart = case_when(
           dt_sum<30 ~ 0,
           (dt_sum>=30)&(dt_sum<60) ~ 1,
           (dt_sum>=60)&(dt_sum<90) ~2,
           dt_sum>90 ~ 2,
            TRUE ~ 99))

sel_MHR %>% group_by(mniQuart) %>%
  summarize(n=n())

rehosp <- masterEXC %>% mutate(rehosp = ifelse(substr(event,1,2) %in% "AI",1,0)) %>%
          filter(rehosp == 1) %>%
            group_by(randomID) %>% distinct(randomID)

masterEXC2 <- masterEXC %>% filter((sel_MHR == 0)&(disTIME %in% c(8)))  %>%
              mutate(eventR = ifelse(((substr(event,1,2) %in% c("AI","NS","CP","CC")))&(substr(event,4,5) %in% c("MI","MG","MO")),
                                     paste(substr(event,1,2),"^","MH",sep=""),event))
  
## segment by first quarter and year of discharge
## masterfq <- master %>% filter(disQTR == 1)

age_des <- masterEXC2 %>% 
                filter(order == 1) %>% 
                group_by(serv_MISS3) %>%
                summarize(mns = mean(AGE_FY16),sdev = sd(AGE_FY16)) %>%
                mutate(mns = round(mns,1), sdev = round(sdev,1)) %>%
                pivot_wider(names_from=c(serv_MISS3),values_from=c(mns,sdev)) %>%
                mutate(variable = "age",value="N/A") %>%
                rename("serv_MISS3_0"="mns_0","serv_MISS3_1"="mns_1",
                       "serv_MISS3_2"="mns_2","freq_0"="sdev_0","freq_1"="sdev_1",
                       "freq_2"="sdev_2")

num_dis <- masterEXC2 %>% 
              filter(order == 1) %>%
              pivot_longer(c(disTIME),"variable","value") %>%
              group_by(variable,value,serv_MISS3) %>%
              summarize(n=n()) %>%
              pivot_wider(names_from=serv_MISS3,values_from=n,names_prefix="serv_MISS3_") 
num_dis <- num_dis %>%
            group_by(variable) %>%
              mutate(freq_0 = round(serv_MISS3_0/sum(serv_MISS3_0)*100,1),
                     freq_1 = round(serv_MISS3_1/sum(serv_MISS3_1)*100,1),
                     freq_2 = round(serv_MISS3_2/sum(serv_MISS3_2)*100,1)) 
                     
num_dis$value <- as.character(num_dis$value)


des_cat <- masterEXC2 %>% 
            filter(order == 1) %>%
          pivot_longer(c(Race,Ethnicity,raceEth),"variable","value") %>%
           group_by(variable,value,serv_MISS3) %>%
            summarize(n = n()) %>%
  pivot_wider(names_from=serv_MISS3,values_from=n,names_prefix="serv_MISS3_")
des_cat <- des_cat %>%
      group_by(variable) %>%
       mutate(freq_0 = round(serv_MISS3_0/sum(serv_MISS3_0)*100,1),
         freq_1 = round(serv_MISS3_1/sum(serv_MISS3_1)*100,1),
         freq_2 = round(serv_MISS3_2/sum(serv_MISS3_2)*100,1)
      )
            
## append the data
desc <- age_des %>%
  bind_rows(.,num_dis) %>%
  bind_rows(.,des_cat) ##%>%
##  bind_rows(.,des_full)

desc <- desc %>%
    mutate(ctsPct_0 = paste(serv_MISS3_0," (",freq_0,")",sep=""),
           ctsPct_1 = paste(serv_MISS3_1," (",freq_1,")",sep=""),
           ctsPct_2 = paste(serv_MISS3_2," (",freq_2,")",sep="")) %>%
  select(variable,value,ctsPct_0,ctsPct_1,ctsPct_2)

## description by events
des_full <- masterEXC2 %>% 
            filter(dt_sum <= 7) %>%
            select(randomID,eventR,serv_MISS3) %>%
            group_by(eventR,serv_MISS3) %>%
            summarize(n=n()) %>%
  pivot_wider(names_from=serv_MISS3,values_from=n,names_prefix="serv_MISS3_")
des_full <- des_full %>% mutate(serv_MISS3_0 = ifelse(is.na(serv_MISS3_0),0,serv_MISS3_0),
                                serv_MISS3_1 = ifelse(is.na(serv_MISS3_1),0,serv_MISS3_1),
                                serv_MISS3_2 = ifelse(is.na(serv_MISS3_2),0,serv_MISS3_2))

des_full <- des_full %>%
          group_by() %>%
          mutate(freq_0 = round(serv_MISS3_0/sum(serv_MISS3_0)*100,1),
                freq_1 = round(serv_MISS3_1/sum(serv_MISS3_1)*100,1),
                freq_2 = round(serv_MISS3_2/sum(serv_MISS3_2)*100,1)) %>%
        mutate(ctsPct_0 = paste(serv_MISS3_0," (",freq_0,")",sep=""),
         ctsPct_1 = paste(serv_MISS3_1," (",freq_1,")",sep=""),
         ctsPct_2 = paste(serv_MISS3_2," (",freq_2,")",sep="")) %>%
  select(eventR,ctsPct_0,ctsPct_1,ctsPct_2)

des_fullEV <- masterEXC2 %>% ##filter(days90 == 1) %>%
## mutate(mhenv = ifelse(substr(event,4,5) %in% c("MH","MI","MG","MO"),1,0)) %>%
  select(randomID,event,serv_MISS3,dt_ID,order) %>%
  group_by(randomID,serv_MISS3) %>%
  summarize(mns_dt=mean(dt_ID),max_ord=max(order))
des_fullEV2 <- des_fullEV %>%
  group_by(serv_MISS3) %>%
  summarize(mns_dt2=mean(mns_dt),sdev_dt=sd(mns_dt),mns_ord=mean(max_ord),sdev_ord=sd(max_ord)) %>%
  mutate(sum_dt=paste(round(mns_dt2,1)," (",round(sdev_dt,1),")",sep=""),
         ord_dt=paste(round(mns_ord,1)," (",round(sdev_ord,1),")",sep="")) %>%
  select(serv_MISS3,sum_dt,ord_dt)

des_fullEV <- masterEXC2 %>% ##filter(days90 == 1) %>%
  mutate(mhenv = ifelse(substr(event,4,5) %in% c("MH","MI","MG","MO"),1,0)) %>%
  select(randomID,event,serv_MISS3,dt_ID,order,mhenv) %>%
  group_by(randomID,serv_MISS3,mhenv) %>%
  summarize(mns_dt=mean(dt_ID),max_ord=max(order))
des_fullEV2 <- des_fullEV %>%
  group_by(serv_MISS3,mhenv) %>%
  summarize(mns_dt2=mean(mns_dt),sdev_dt=sd(mns_dt),mns_ord=mean(max_ord),sdev_ord=sd(max_ord)) %>%
  mutate(sum_dt=paste(round(mns_dt2,1)," (",round(sdev_dt,1),")",sep=""),
         ord_dt=paste(round(mns_ord,1)," (",round(sdev_ord,1),")",sep="")) %>%
  select(serv_MISS3,sum_dt,ord_dt,mhenv)
##pivot_wider(names_from=serv_MISS3,values_from=n,names_prefix="serv_MISS3_")

p <- ggplot(masterEXC2,aes(x=dt_ID)) + geom_bar() +
  ggtitle("Distribution of Time between Events") +
  labs(x="Time between Events",y="Counts") 
 ##xlim(100,3500)
p

p <- ggplot(des_fullEV,aes(x=max_ord)) + geom_bar() +
  ggtitle("Distribution of Number of Events per Person") +
  labs(x="Events per Person",y="Count") 
##xlim(100,3500)
p

## Wilcoxon rank sum test used to test for notable differences between medians
## signed rank test used to test for notable differences between before and after the intervention
## mean and standard deviation presented, t-test used
## fischer's exact test used and to form a cell-count less than 5
## chi-square test for categorical variables
## age
temp <- masterEXC2 %>% 
        filter((days90 == 1)&(serv_MISS3 != 1)&(order == 1))
tbl <-  table(temp$raceEth,temp$serv_MISS3)
fisher.test(tbl)
temp2 <- masterEXC2 %>% 
        filter((days90 == 1)&(serv_MISS3 != 0)&(order == 1))
tbl2  <-  table(temp2$raceEth,temp2$serv_MISS3)
fisher.test(tbl2)
A <- masterEXC2 %>% filter((days90 == 1)&(serv_MISS3 == 0)&(order == 1)) %>% select(AGE_FY16)
B <- masterEXC2 %>% filter((days90 == 1)&(serv_MISS3 == 1)&(order == 1)) %>% select(AGE_FY16)
C <- masterEXC2 %>% filter((days90 == 1)&(serv_MISS3 == 2)&(order == 1)) %>% select(AGE_FY16)
t.test(A$AGE_FY16,C$AGE_FY16)
t.test(B$AGE_FY16,C$AGE_FY16)


#############################################################################
## cluster the sequences

first_date = ymd(20000101)
# 
# df_ncfw <-
#   master %>%
#  select(randomID,dt_sum,event,order,subs,serv_MISS) %>%
#   filter((subs == 1)&(serv_MISS == 1)) %>%
#   arrange(randomID,order) %>%
#   group_by(randomID) %>%
#   mutate(period = (first_date + dt_sum) %>% as.character()) %>%
#   select(id = randomID,
#          period,
#          event)
# 
# ## This runs the week clustering.
# df_ncfw %>%
#   aggregate_sequences(format = "%Y-%m-%d",unit="week",n_units=1) %>%
#   cluster_knn(k=5) %>%
#   filter_pattern(threshold = 0.3, pattern_name = "variation") %>% 
#   filter_pattern(threshold = 0.5, pattern_name = "consensus") %>% 
#   filter_pattern(threshold = 0.7, pattern_name = "consprioritty") %>% 
#   generate_reports(output_directory = "./output/")
#   
# ## This runs the month clustering.
# df_ncfw %>%
#   aggregate_sequences(format = "%Y-%m-%d",unit="month",n_units=1) %>%
# cluster_knn(k=5) %>%
#   filter_pattern(threshold = 0.3, pattern_name = "variation") %>% 
#  filter_pattern(threshold = 0.5, pattern_name = "consensus") %>% 
#  filter_pattern(threshold = 0.7, pattern_name = "consprioritty") %>% 
#   generate_reports(output_directory = "./output/")

################################################################
## run the same for first 90 days only with clustering by week


df_ncfw_90d <-
  masterEXC2 %>%
  select(randomID,dt_sum.x,eventR,order,serv_MISS,days90) %>%
  filter(serv_MISS == 1) %>%
  arrange(randomID,order) %>%
  group_by(randomID) %>%
  mutate(period = (first_date + dt_sum.x) %>% as.character()) %>%
  select(id = randomID,
         period,
         event = eventR)

## Cluster by week.
temp <- df_ncfw_90d %>%
  aggregate_sequences(format = "%Y-%m-%d",unit="week",n_units=1) %>%
  cluster_knn(k=3) %>%
  filter_pattern(threshold = 0.3, pattern_name = "variation") %>% 
  filter_pattern(threshold = 0.5, pattern_name = "consensus")## %>% 
 ##filter_pattern(threshold = 0.7, pattern_name = "consprioritty")## %>% 
 ## generate_reports(output_directory = "./output/")

temp2 <-  temp %>%
  mutate_all(.tbl=.,funs(replace(.,lengths(.)==0,list(list(elements=,element_weights=,)))))
temp2 %>% generate_reports(output_directory = "./output/")

## cluster 2: 3627, 6377, 9859, 10146, 2975, 2554, 3402, 5883, 2641, 
## 5986, 11618, 11682, 8182, 10301, 3704, 8857, 11453, 1256, 3793, 
## 3567, 8435, 8714, 1109, 8714, 1109, 5162, 7663, 8159, 9402, 6455,
## 3146, 6220, 7774 2837, 591, 5164, 3707, 4398, 11010, 7654, 12022,
## 6292, 9334, 8348, 1584, 6109, 6610, 2088
clust2 <- masterEXC2 %>% filter(randomID %in% c(3627, 6377, 9859, 10146, 2975, 2554, 3402, 5883, 2641, 
                                                5986, 11618, 11682, 8182, 10301, 3704, 8857, 11453, 1256, 3793, 
                                                567, 8435, 8714, 1109, 8714, 1109, 5162, 7663, 8159, 9402, 6455,
                                                3146, 6220, 7774, 2837, 591, 5164, 3707, 4398, 11010, 7654, 12022,
                                                6292, 9334, 8348, 1584, 6109, 6610, 2088)) %>%
  select(randomID,eventR,dt_sum) %>%
  mutate(dt_week = ceiling((dt_sum+1)/7))

## cluster 2: 3722, 1543, 3104, 561, 4853, 5372, 6785, 4388, 10295,
## 1157, 3491, 5868, 8474, 7387, 9148, 8525, 5819, 8933, 2615,
## 2938, 3760, 6314, 7386, 11251, 3621, 3030, 4334, 4308, 155

clust2 <- masterfq %>% filter(days90 ==1, randomID %in% c(3722, 1543,3104,561,
                                              4853,5372,6875,4388,10295,
                                              1157,3491,5868,8474,7387,
                                              9148,8525,5819,8933,2615,
                                              2938,3760,6314,7386,11251,
                                              3621,3030,4334,4308,155)) %>%
                      select(randomID,event,dt_sum)

################################################################
## rerun by race/ethnicity:
##  "NOT HISPANIC WHITE", "HISPANIC OR LATINO", "NOT HISPANIC AFRICAN AMERICAN"

############################
## non hispanic white 

df_ncfw_nhw <-
  master %>%
  select(randomID,dt_sum,event,order,serv_MISS,raceEth,days90) %>%
  filter((serv_MISS == 1)&(raceEth == "NOT HISPANIC WHITE")&(days90 == 1)) %>%
  arrange(randomID,order) %>%
  group_by(randomID) %>%
  mutate(period = (first_date + dt_sum) %>% as.character()) %>%
  select(id = randomID,
         period,
         event)

## This runs the set for week.
temp <- df_ncfw_nhw %>%
  aggregate_sequences(format = "%Y-%m-%d",unit="week",n_units=1) %>%
  cluster_knn(k=5) %>%
  filter_pattern(threshold = 0.3, pattern_name = "variation") %>% 
  filter_pattern(threshold = 0.5, pattern_name = "consensus") ##%>% 
##  filter_pattern(threshold = 0.7, pattern_name = "consprioritty") %>% 
##  generate_reports(output_directory = "./output/")

temp %>% save_alignment() %>% write_file("alignments.csv")

############################
## hispanic or latino 

df_ncfw_lat <-
  master %>%
  select(randomID,dt_sum,event,order,serv_MISS,raceEth,days90) %>%
  filter((serv_MISS == 1)&(raceEth == "HISPANIC OR LATINO")&(days90 == 1)) %>%
  arrange(randomID,order) %>%
  group_by(randomID) %>%
  mutate(period = (first_date + dt_sum) %>% as.character()) %>%
  select(id = randomID,
         period,
         event)

## This runs the set for week.
df_ncfw_lat %>%
  aggregate_sequences(format = "%Y-%m-%d",unit="week",n_units=1) %>%
  cluster_knn(k=5) %>%
  filter_pattern(threshold = 0.3, pattern_name = "variation") %>% 
  filter_pattern(threshold = 0.5, pattern_name = "consensus") %>% 
##  filter_pattern(threshold = 0.7, pattern_name = "consprioritty") %>% 
  generate_reports(output_directory = "./output/")

############################
## non hispanic african american
df_ncfw_nhafam <-
  master %>%
  select(randomID,dt_sum,event,order,serv_MISS,raceEth,days90) %>%
  filter((serv_MISS == 1)&(raceEth == "NOT HISPANIC AFRICAN AMERICAN")&(days90 == 1)) %>%
  arrange(randomID,order) %>%
  group_by(randomID) %>%
  mutate(period = (first_date + dt_sum) %>% as.character()) %>%
  select(id = randomID,
         period,
         event)

## This runs the set for week.
df_ncfw_nhafam %>%
  aggregate_sequences(format = "%Y-%m-%d",unit="week",n_units=1) %>%
  cluster_knn(k=5) %>%
  filter_pattern(threshold = 0.3, pattern_name = "variation") %>% 
  filter_pattern(threshold = 0.5, pattern_name = "consensus") %>% 
##  filter_pattern(threshold = 0.7, pattern_name = "consprioritty") %>% 
  generate_reports(output_directory = "./output/")

#######################################################################
## follow-up on hypotheses
masterEXC3 <- rename(masterEXC2, dt_sum = dt_sum.x)
df_hypoth <-
  masterEXC3 %>% 
  select(randomID,dt_sum,eventR,order,serv_MISS,days90,dt_ID) %>%
  filter(serv_MISS == 1)

## develop counts for hypothesis one
count_MH <- df_hypoth %>% filter(eventR == c("AT^XM")) %>% 
              group_by(randomID)
length(unique(count_MH$randomID))

temp <- count_MH %>% filter(dt_sum < 8)
length(unique(temp$randomID))

temp <- count_MH %>% group_by(randomID) %>%
            summarize(n=n(),msnss=mean(dt_ID))
temp2 <- temp %>% summarize(mns=mean(n),sdev=sd(n),mnsdt=mean(msnss),sddt=sd(msnss))

## average time re-attend after first week
temp <- count_MH %>% filter(dt_sum > 7) %>% group_by(randomID) %>%
        summarize(n=n(),mnss=mean(dt_ID))
temp2 <- temp %>% summarize(mns=mean(n),sdv=sd(n),mnsdt=mean(mnss),sddt=sd(mnss))

## average time re-attend after missing a mental health appt in first seven days
temp <- count_MH %>% 
        filter(dt_sum > 7) %>%
        group_by(randomID) %>%
        mutate(dt_sum7 = dt_sum - 7) %>%
         slice(which.min(order)) %>% 
        as.data.frame(.) %>% ungroup() %>%
        summarize(mns=mean(dt_sum7),sdev=sd(dt_sum7))

## cluster 2 those who re-rengaged with mental health treatment in 90 day period
count_reMH <- df_hypoth %>% filter(eventR %in% c("AT^MI","AT^MG","AT^MO"))
length(unique(count_reMH$randomID))

temp <- count_reMH %>% group_by(randomID) %>% 
        summarize(n=n(),mnss=mean(dt_ID)) %>% ungroup() %>% 
      summarize(mns=mean(n),std=sd(n),mnsID=mean(mnss),sdID=sd(mnss))

## cluster 3 those who had minimal re-engagement or were lost to follow-up
## create variable for: number of follow-up mental health visits, <1,1-4,4-8,>8
temp <- count_reMH %>% group_by(randomID) %>% 
        summarize(n=n())
tempct <- dplyr::right_join(temp,df_hypoth,by=c("randomID"))
tempct2 <- tempct %>% mutate(n = ifelse(is.na(n),0,n)) %>%
          mutate(vsts = case_when(
            (n == 0) ~ 0,
            (n == 1) ~ 1,
            (n > 1)&(n <= 4) ~ 2,
            (n > 4)&(n <= 8) ~ 3,
            (n > 8) ~ 4,
            TRUE ~ 99)) %>% as.data.frame()
tempct3 <- data.frame(randomID=tempct2$randomID,vsts=tempct2$vsts,
                         order=tempct2$order,n=tempct2$n)
##tempct3 <- as.data.frame(tempct2)
temp <- tempct2 %>% group_by(randomID) %>%
              slice(which.min(order)) %>% ungroup() %>%
              group_by(vsts) %>% summarize(nt=n())


  