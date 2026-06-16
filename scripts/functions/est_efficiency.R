#' @title est_efficiency
#'
#' @description Estimate efficiency for traps using mark-recapture efficiency trials. 
#' Called by est_passage function.
#' Takes the following steps:
#' 1) Align release events with trap visits within 36 hours post-release.
#' 2) Calculate efficiency as recaptures / releases per batch date.
#' 3) Expand to full survey period (all days between min/max trial dates).
#' 4) Calculate mean recapture time weighted by number of recaps.
#' 5) Pass trial efficiency data to model_efficiency().
#' 
#' @param release_data Data set with records of fish marked and released for efficiency trials.
#' 
#' @param recapture_data Data set with records of marked fish recaptured for efficiency trials.
#' 
#' @param impute_all User input of dictating if all efficiency estimates should be imputed. Defaults to FALSE. If TRUE, will replace estimates made with actual
#' efficiency trial data with imputed values.
#' 
#' @param visit_data Input trap visit data frame.
#' 
#' @param survey_start Input survey season start date.
#' 
#' @param survey_end Input survey season end date.
#' 
#' @param min_sample_size Minimum number of efficiency trials required in order to fit binomial model to impute efficiency. Defaults to 10.
#' 
#' @return eff_modeled List containing efficiency model fit results and data with imputed efficiency values.

est_efficiency<-function(release_data,
                         recapture_data,
                         visit_data,
                         impute_all=F,
                         survey_start,
                         survey_end,
                         min_sample_size=10){
  
  #rename
  visits<-visit_data
  releases<-release_data
  recaps<-recapture_data
  
  #fix dates
  visits<-fix_dates(visits,"visit_time")
  releases<-fix_dates(releases,"release_time")
  recaps<-fix_dates(recaps,"visit_time")
  
  #get sampling time and trap_ID_decimal for visits
  visits<-calculate_sampling_time(visits)
  
  #change datetime names to differentiate visit vs release times
  releases<-releases%>%
    rename(release_datetime=date_time)
  visits<-visits%>%
    rename(visit_datetime=date_time)
  recaps<-recaps%>%
    rename(visit_datetime=date_time)
  
  #add the trap_ID value to have a unique ID for each subsite x river
  #this is in case we have data for multiple sites
  recaps$trap_ID<-paste(recaps$site_name,
                                recaps$subsite_name)
  
  #filter by those records marked for analysis
  release_use <- releases %>%
    filter(include_analysis=="Yes")
  
  #create new dataframe of all visits within default 36 hours of a release
  #for each release, find the relevant visits using the 
  #"find_visits_near_release" function
  visits_within_36_hours <- release_use %>%
    mutate(visits_in_window = map(release_datetime, 
                                  find_visits_near_release, 
                                  visit_df = visits)) %>%
    unnest(visits_in_window) %>%
    select(release_ID,release_datetime, trap_ID,trap_ID_decimal,
           trap_visit_ID, visit_datetime,n_released) %>%
    mutate(hours_after_release = as.numeric(difftime(visit_datetime,
                                                     release_datetime,
                                                     units = "hours")))
  #determine number of hours between release and first and last visits within
  #36 hours of release (usually just one visit)
  visits_within_36_hours<-visits_within_36_hours%>%
    group_by(release_ID, trap_ID,trap_ID_decimal,release_datetime,n_released) %>%
    summarize(HrsToFirstVisitAfter = as.numeric(difftime(min(visit_datetime),
                                                         first(release_datetime),
                                                         units = "hours")),
              HrsToLastVisitAfter = as.numeric(difftime(max(visit_datetime),
                                                        last(release_datetime),
                                                        units = "hours")))
  
  #join visits to recap data
  trial_data<-visits_within_36_hours%>%
    left_join(select(recaps,
                     trap_visit_ID,release_ID,n,
                     trap_ID,visit_datetime))
  
  #calculate mean recapture time weighted by number of recaps
  trial_data <- trial_data %>%
    group_by(release_ID, trap_ID,trap_ID_decimal,
             HrsToFirstVisitAfter,HrsToLastVisitAfter) %>%
    summarize(
      n_released = first(n_released),
      release_datetime = first(release_datetime),
      n_recaps = sum(n, na.rm = TRUE),
      meanRecapTime = sum(as.numeric(visit_datetime) * n) / sum(n),
      .groups = 'drop'
    )
  
  #set meanRecapTime for trials with no recaptures
  trial_data<-trial_data%>%
    mutate(meanRecapTime=ifelse(n_recaps==0,
                                as.numeric(release_datetime)+((HrsToFirstVisitAfter+HrsToLastVisitAfter)/2)*3600,
                                meanRecapTime))
  
  #set time/date format for meanRecapTime
  trial_data$meanRecapTime<-as.POSIXct(trial_data$meanRecapTime,
                                       origin="1970-01-01",
                                       tz="America/Los_Angeles")
  
  #assign batch_date to trial_data
  #this will allow us to match catch data to trial batch_dates
  trial_data<-assign_batch_date(df=trial_data,
                              time_field="meanRecapTime")
  
  #sometimes (2018-10-25) we have releases on same day but with different markers
  #so here we combine
  trial_data<-trial_data%>%
    group_by(trap_ID,trap_ID_decimal,meanRecapTime,batch_date)%>%
    summarise(n_released=sum(n_released),
              n_recaps=sum(n_recaps))
  
  #estimate efficiency for each releaseID trap_ID combo
  trial_data$n_released[ trial_data$n_released <= 0] <- NA 
  trial_data$efficiency <- (trial_data$n_recaps)/(trial_data$n_released)
  trial_data <- trial_data[ !is.na(trial_data$efficiency), ]
  
  #set survey date range based on user input
  start_date<-as.POSIXct(survey_start,
                         origin="1970-01-01",
                         tz="America/Los_Angeles")
  end_date<-as.POSIXct(survey_end,
                       origin="1970-01-01",
                       tz="America/Los_Angeles")
  
  trial_data<-trial_data%>%
    filter(batch_date<=end_date & batch_date>=start_date)
  
  season<-seq(start_date,end_date,by="days")
  
  #figure out which days have efficiency data
  survey_period <- expand.grid(trap_ID_decimal=sort(unique(trial_data$trap_ID_decimal)),
                               batch_date=format(season,"%Y-%m-%d"),
                               stringsAsFactors=F)
  survey_period$batch_date<-as.POSIXct(survey_period$batch_date,
                                       origin="1970-01-01",
                                       tz="America/Los_Angeles")
  eff<-survey_period%>%
    left_join(trial_data,by=c("trap_ID_decimal","batch_date"))
  
  #lets save eff inputs
  eff$unimputed_efficiency<-eff$efficiency
  #next build eff_model and run here
  
  eff_modeled<-model_efficiency(efficiency_data=eff,
                        impute_all=impute_all,
                        min_sample_size=min_sample_size)
  
  eff_modeled$results<-unique(eff_modeled$results)%>%
    left_join(dplyr::select(eff,
                     batch_date,trap_ID_decimal,
                     unimputed_efficiency))
  
  return(eff_modeled)
}
