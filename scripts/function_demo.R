#Function demo script
#for testing the functions
rm( list = ls()) #clear env

library(tidyverse)

###############################
#example data from feather river
###############################

catch<-read.csv("data/test_data/demo_catch.csv")
release<-read.csv("data/test_data/demo_release.csv")
recapture<-read.csv("data/test_data/demo_recapture.csv")
visits<-read.csv("data/test_data/demo_visit.csv")

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

#testing grounds
survey_start<-"2022-01-19"
survey_end<-"2022-06-22"
target_species<-"Chinook salmon"
target_run<-NA
target_run<-"Fall"

summarize.by="week"
impute_all=F
bootstrap=T
file.name="test"

test_pass<-est_passage(catch,visits,
                       release,recapture,
                       summarize.by="week", 
                       impute_all=F,
                       bootstrap=T,
                       survey_start,survey_end,
                       target_species,
                       target_run,
                       file.name="test")
plot(test_pass$p_passage)

