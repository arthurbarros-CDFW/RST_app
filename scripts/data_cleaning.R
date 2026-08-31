#data cleaning scripts
#script to clean raw rst data files from diffrent sites
rm( list = ls()) #clear env

library(tidyverse)
library(readxl)
library(scales)

###############################
#battle creek
###############################

recapture<-read.csv("data/battle/recapture.csv")
release<-read.csv("data/battle/release.csv")

#visit data
visit<-read.csv("data/battle/visit.csv")

visit<-visit%>%
  rename(site_name=station_code,
         trap_visit_id=sample_id)
visit<-visit%>%
  mutate(subsite_name=ifelse(is.na(ubc_site),1,ubc_site),
         visit_time=as.POSIXct(paste(sample_date, sample_time),
                               format="%Y-%m-%d %H:%M:%S"))

#catch data
catch<-read.csv("data/battle/catch.csv")

catch<-catch%>%left_join(visit)
catch<-catch%>%
  rename(at_capture_run=fws_run,
         )