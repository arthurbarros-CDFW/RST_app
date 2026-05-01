#Function demo script
#for testing the functions
rm( list = ls()) #clear env

library(tidyverse)
library(readxl)
library(scales)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

###############################
#BUTTE
###############################

catch<-read_xlsx("data/butte/butte_catch.xlsx")
release<-read_xlsx("data/butte/butte_release.xlsx")
recapture<-read_xlsx("data/butte/butte_recapture.xlsx")
visits<-read_xlsx("data/butte/butte_trap.xlsx")

#testing grounds
survey_start<-"2021-01-19"
survey_end<-"2021-06-22"
target_species<-"Chinook salmon"
#target_run<-"Fall" turn off because run not recorded properly

butte_pass<-est_passage(catch,
                        visits,
                        release,
                        recapture,
                       summarize.by="week", 
                       impute_all=F,
                       bootstrap=T,
                       survey_start,survey_end,
                       target_species,
                       file.name="butte")


###############################
#KNIGHTS LANDING
###############################


catch<-read.csv("data/knights landing/knights_landing_catch.csv")
release<-read.csv("data/knights landing/knights_landing_release.csv")
recapture<-read.csv("data/knights landing/knights_landing_recapture.csv")
visits<-read.csv("data/knights landing/knights_landing_trap.csv")

#testing grounds
survey_start<-"2018-01-01"
survey_end<-"2021-06-22"
target_species<-"Chinook salmon"
target_run<-"Fall" 
summarize.by="week"
impute_all=F
bootstrap=T
file.name="knights"

knights_pass<-est_passage(catch,
                        visits,
                        release,
                        recapture,
                        summarize.by="week", 
                        impute_all=F,
                        bootstrap=T,
                        survey_start,survey_end,
                        target_species,
                        file.name="knights")


