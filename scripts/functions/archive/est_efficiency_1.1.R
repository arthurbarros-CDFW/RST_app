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
                         min_sample_size=10,
                         use_discharge=FALSE){
  
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
  
  if(nrow(visits_within_36_hours) == 0) {
    stop("No visits found within 36 hours of releases")
  }
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
  
  if(nrow(trial_data) == 0) {
    stop("No efficiency trials could be constructed from the data")
  }
  
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
  
  #sometimes we have releases on same day but with different markers
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
  
  if(nrow(trial_data) == 0) {
    stop("No efficiency trials found within the specified survey dates")
  }
  
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
  
  if(use_discharge == TRUE){
    visits <- assign_batch_date(df = visits,
                                time_field = "visit_datetime")
    
    #aggregate discharge to daily values per trap
    env_covars <- visits %>%
      ungroup() %>%
      group_by(trap_ID_decimal, batch_date) %>%
      summarise(discharge = mean(discharge, na.rm = TRUE), .groups = "drop") %>%
      select(trap_ID_decimal, batch_date, discharge)
    
    eff <- eff %>%
      left_join(env_covars, by = c("trap_ID_decimal", "batch_date"))
    
    eff$discharge_scaled <- as.numeric(scale(eff$discharge))
  }
  
  eff_modeled<-model_efficiency(efficiency_data=eff,
                                impute_all=impute_all,
                                min_sample_size=min_sample_size,
                                use_discharge = use_discharge)
  
  eff_modeled$results<-unique(eff_modeled$results)%>%
    left_join(dplyr::select(eff,
                            batch_date,trap_ID_decimal,
                            unimputed_efficiency))

  #filter out NA discharge
  if(use_discharge){
    discharge_data <- eff_modeled$results %>%
      filter(!is.na(discharge) & !is.nan(discharge))
    scale_factor <- max(eff_modeled$results$efficiency, na.rm = TRUE) / 
      max(eff_modeled$results$discharge, na.rm = TRUE)
  }
  
  p_eff <- ggplot()
  
  #add conditional layers
  if (impute_all == FALSE) {
    p_eff <- p_eff +
      geom_point(data = eff_modeled$results, 
                 aes(x = batch_date, y = efficiency_imputed, 
                     color = "Imputed Efficiency",
                     shape = "Imputed Efficiency")) +
      geom_point(data = eff_modeled$results, 
                 aes(x = batch_date, y = unimputed_efficiency, 
                     color = "Unimputed Efficiency",
                     shape = "Unimputed Efficiency"))
  } else if (impute_all == TRUE) {
    p_eff <- p_eff +
      geom_point(data = eff_modeled$results, 
                 aes(x = batch_date, y = efficiency_imputed, 
                     color = "Imputed Efficiency",
                     shape = "Imputed Efficiency"))
  }
  
  #add discharge points
  if (use_discharge == TRUE) {
    p_eff <- p_eff +
      geom_point(data = discharge_data, 
                 aes(x = batch_date, y = discharge * scale_factor, 
                     color = "Discharge",
                     shape = "Discharge"), 
                 alpha = 0.6, size = 1.5) +
      scale_color_manual(
        name="Values",
        values = c("Imputed Efficiency" = "#00BFC4", 
                   "Unimputed Efficiency" = "#F8766D",
                   "Discharge" = "black")
      ) +
      scale_shape_manual(
        name="Values",
        values = c("Imputed Efficiency" = 15,    
                   "Unimputed Efficiency" = 16,  
                   "Discharge" = 21)
      ) +
      scale_y_continuous(
        name = "Efficiency",
        sec.axis = sec_axis(
          ~ . / scale_factor,
          name = "Discharge (cfs)"
        )
      )
  } else{
    p_eff <- p_eff +
      scale_color_manual(
        name="Values",
        values = c("Imputed Efficiency" = "#00BFC4", 
                   "Unimputed Efficiency" = "#F8766D")
      ) +
      scale_shape_manual(
        name="Values",
        values = c("Imputed Efficiency" = 15,    
                   "Unimputed Efficiency" = 16)
      ) +
      scale_y_continuous(
        name = "Efficiency"
      )
  }
  
  p_eff <- p_eff +
    theme_bw() +
    facet_wrap(. ~ trap_ID_decimal)
  
  eff_modeled$p_eff<-p_eff
  return(eff_modeled)
}
