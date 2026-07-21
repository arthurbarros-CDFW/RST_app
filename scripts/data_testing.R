#Function demo script
#for testing the functions
rm( list = ls()) #clear env

library(tidyverse)
library(readxl)
library(scales)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

###############################
#Efficiency ~ discharge testing
###############################

catch<-read.csv("data/test_data/demo_catch.csv")
releases<-read.csv("data/test_data/demo_release.csv")
recaps<-read.csv("data/test_data/demo_recapture.csv")
visits<-read.csv("data/test_data/demo_visit.csv")

#testing grounds
survey_start<-"2022-01-19"
survey_end<-"2022-06-22"
target_species<-"Chinook salmon"
#target_run<-"Fall" turn off because run not recorded properly

impute_all=F
min_sample_size=10

test_results<-est_efficiency(release_data=releases,
                                  recapture_data=recaps,
                                  visit_data=visits,
                                  impute_all=F,
                                  survey_start,
                                  survey_end,
                                  min_sample_size=10,
                                  use_discharge=TRUE)

###############################
#Passage ~ discharge testing
###############################

catch<-read.csv("data/test_data/demo_catch.csv")
release<-read.csv("data/test_data/demo_release.csv")
recapture<-read.csv("data/test_data/demo_recapture.csv")
visits<-read.csv("data/test_data/demo_visit.csv")

#testing grounds
survey_start<-"2022-01-19"
survey_end<-"2022-06-22"
target_species<-"Chinook salmon"
target_run<-"Fall"

summarize.by="week"
impute_all=F
bootstrap=T
file.name="test"
use_discharge=TRUE

test_pass<-est_passage(catch,visits,
                       release,recapture,
                       summarize.by="week", 
                       impute_all=F,
                       bootstrap=T,
                       survey_start,survey_end,
                       target_species,
                       target_run,
                       file.name="test",
                       use_discharge = use_discharge)
plot(test_pass$p_eff)

