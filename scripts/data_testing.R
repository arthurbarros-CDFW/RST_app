#Function demo script
#for testing the functions
rm( list = ls()) #clear env

library(tidyverse)
library(readxl)
library(scales)

#Load functions
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

##############################
#dealing with visit_type and fish_processed
##############################
visits<-read.csv("data/test_data/demo_visit.csv")
v<-visits
#fix dates
v<-fix_dates(v,"visit_time")

#set batch dates
v<-assign_batch_date(v,"date_time")


#get sampling time
v<-calculate_sampling_time(v)

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
target_run<-"Fall"

impute_all=F
min_sample_size=10
release_data=releases
recapture_data=recaps
visit_data=visits
catch_data = catch
use_discharge=TRUE
max.df.spline=4
sum.by="week"

test_catch_results<-est_catch(target_species = target_species,
                              target_run = target_run,
                              catch_data = catch,
                              visit_data = visits,
                              survey_start = survey_start,
                              survey_end = survey_end)

test_eff_results<-compare_efficiency_models(release_data=releases,
                                  recapture_data=recaps,
                                  visit_data=visits,
                                  impute_all=F,
                                  survey_start,
                                  survey_end,
                                  min_sample_size=10,
                                  use_discharge=use_discharge)


eff.data.unimputed<-test_eff_results$eff
test_eff_types<-test_eff_results$eff.type

test_selected_models<-list(
  "Lower Feather River RST RL"=test_eff_results$candidate_models$`Lower Feather River RST RL`$`ROM+1 (constant efficiency)`,
  "Lower Feather River RST RR"=test_eff_results$candidate_models$`Lower Feather River RST RR`$temporal_df4
)

#in app we pull this from the select models reactive chunk, it doesn't need to 
#pass to imputation, just passage estimate for bootstrapping
test_selected_types<-list(
  "Lower Feather River RST RL"=1,
  "Lower Feather River RST RR"=5
)

model_names_test<-list(
  "Lower Feather River RST RL"="ROM+1 (constant efficiency)",
  "Lower Feather River RST RR"="discharge_temporal_df4"
)
  

test_eff<-impute_efficiency(efficiency_data=eff.data.unimputed,
                            selected_models=test_selected_models,
                            model_names=model_names_test,
                            impute_all=impute_all)

catch.results = test_catch_results
eff.results = test_eff
summarize.by = sum.by
catch.fits = test_catch_results$models
eff.fits = test_eff$models
eff.types=test_selected_types
bootstrap=T

test_passage<-est_passage(catch.results=catch.results,
                          eff.results=eff.results,
                          summarize.by=summarize.by,
                          survey_start=survey_start,
                          survey_end=survey_end,
                          catch.fits=catch.fits,
                          eff.fits=eff.fits,
                          eff.types=eff.types,
                          target_species=target_species,
                          target_run=target_run,
                          min_sample_size = 10)

