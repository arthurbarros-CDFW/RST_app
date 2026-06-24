#' @title est_catch
#'
#' @description Estimate daily catch per trap and impute missing catch values. Called by est_passage function.
#' Takes the following steps:
#' 1) Join catch data to visit data.
#' 2) Group and sum catch by trap_ID and batch_date.
#' 3) Expand date range to include all days between first and last sampling.
#' 4) Create “empty” records for not fishing periods.
#' 5) Calls model_catch() for imputation of missing catch values.
#' 
#' @param target_species Input of target species, defaults to "Chinook Salmon".
#' 
#' @param target_run Input of target run.
#' 
#' @param catch_data Input catch data frame.
#' 
#' @param visit_data Input trap visit data frame.
#' 
#' @param survey_start Input survey season start date.
#' 
#' @param survey_end Input survey season end date.
#' 
#' @return final_catch List containing catch model fit results and catch data with imputed catch values.

est_catch<-function(target_species,
                    target_run=NA,
                    catch_data,
                    visit_data,
                    survey_start,
                    survey_end){
  v<-visit_data
  c<-catch_data
  #fix dates
  c<-fix_dates(c,"visit_time")
  v<-fix_dates(v,"visit_time")
  
  #set batch dates
  c<-assign_batch_date(c,"date_time")
  v<-assign_batch_date(v,"date_time")
  
  
  #get sampling time
  v<-calculate_sampling_time(v)
  
  #let's add "trapID" value to catch
  c$trap_ID<-paste(c$site_name,
                            c$subsite_name)
  
  #filter and join catch to visit data
  
  #if no target run just all runs
  if(is.na(target_run)){
    target_catch<-c%>%filter(
      common_name==target_species
    )
  }else{
    target_catch<-c%>%filter(
      common_name==target_species &
        at_capture_run==target_run
    )
  }
  target_catch<-target_catch%>%
    left_join(select(v,trap_visit_ID,
                     trap_ID,
                     include_catch,
                     batch_date,
                     start_time,
                     end_time))
  
  #sum catch by batch date and trap
  est_catch<-target_catch%>%
    group_by(trap_ID,batch_date)%>%
    summarize(total_catch=sum(n))
  
  visit_catch<-v%>%
    left_join(est_catch)%>%
    mutate(
      #if trap was fishing and there is no catch
      #set total catch to 0
      total_catch=ifelse(is.na(total_catch) & trap_status=="fishing",0,total_catch)
    )%>%
    filter(include_catch=="Yes")
  
  batch_date_catch<-visit_catch%>%
    group_by(batch_date,trap_ID,trap_ID_decimal,
             trap_status,start_time,end_time)%>%
    summarise(total_sample_minutes=sum(total_sample_minutes),
              total_catch=sum(total_catch))
  
  #set survey date range based on user input
  start_date<-as.POSIXct(survey_start,
                         origin="1970-01-01",
                         tz="America/Los_Angeles")
  end_date<-as.POSIXct(survey_end,
                       origin="1970-01-01",
                       tz="America/Los_Angeles")
  
  batch_date_catch<-batch_date_catch%>%
    filter(batch_date<=end_date & batch_date>=start_date)
  
  #next we want to fill in blanks for other dates not showing up in visits (ie days between trap ends and starts)
  #needs to be broken up by trap_ID_decimal
  
  final_catch <- batch_date_catch %>%
    group_by(trap_ID,trap_ID_decimal,start_time,end_time) %>%
    complete(batch_date = seq(min(batch_date), max(batch_date), by = "day"),
             fill = list(trap_status = "Not fishing", 
                         total_sample_minutes =0, 
                         total_catch = NA)) %>%
    fill(trap_ID, .direction = "downup") %>%  #fill missing trap_ID values
    ungroup()
  
  #recalculate total_sample_minutes for Not fishing periods
  final_catch<-final_catch%>%
    mutate(total_sample_minutes=ifelse(
      !is.na(start_time),
      as.numeric(difftime(end_time,start_time,units="mins")),
      NA_real_
    ))
  
  final_catch<-model_catch(final_catch)
  
  return(final_catch)
}
