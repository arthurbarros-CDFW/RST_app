#' @title calculate_sampling_time
#'
#' @description Process raw trap visit data to derive continuous sampling periods.
#' 1) Filters for valid visit types.
#' 2) Calculate start/end times and duration between visits.
#' 3) Classify periods as either "fishing" or "Not fishing".
#' 4) Remove short non-fishing periods less than 30 minutes.
#' 5) Split sampling into segments separated by large gaps greater than gap_threshold_days.
#' 6) Create trap_ID_decimal field as a unique ID for trap sampling segments.
#' 
#' @param visit_data Input trap visit data frame.
#' 
#' @return data Outputs data frame with traps broken up by sampling segments.
#' 
calculate_sampling_time<-function(visit_data){
  #identify visits to include
  data=visit_data
  data<-data%>%
    filter((visit_type%in%c("Start trap & begin trapping",
                           "Continue trapping",
                           "Unplanned restart",
                           "End trapping") &
              fish_processed != "No catch data; fish left in live box") |
             (visit_type=="Start trap & begin trapping" &
                fish_processed == "No catch data; fish left in live box")
    )
  
  #identify time_sample_end as the datetime
  #calculate time_sample_start as the prior visits datetime
  
  #set start and end times for sampling
  data <- data %>%
    arrange(subsite_name, date_time)%>%
    group_by(subsite_name)%>%
    mutate(
      order_start_sample = row_number(),
      end_time = date_time,
      start_time = lag(date_time),
      total_sample_minutes = ifelse(
        !is.na(start_time),
        as.numeric(difftime(end_time, start_time, units = "mins")),
        NA_real_
      )
    ) %>%
    ungroup()
  
  
  data$trap_ID=paste(data$site_name,data$subsite_name)
  
  #figure out periods when fishing vs not fishing
  data <- data %>%
    arrange(subsite_name, date_time)%>%
    group_by(subsite_name)%>%
    mutate(
      trap_status=ifelse(
        lag(visit_type)=="End trapping" | visit_type=="Start trap & begin trapping",
        "Not fishing","fishing"
      )
    )
  
  #filter out "not fishing" periods less than 30 min
  data_out<-data%>%
    filter(
      trap_status=="Not fishing" & total_sample_minutes<=30
    )
  
  data<-data%>%
    filter(!trap_visit_ID %in% unique(data_out$trap_visit_ID))
  
  #next separate sampling groups by long gaps 
  #(default= 7 days, 10080 minutes)
  gap_threshold_minutes=gap_threshold_days*24*60
  
  data<-data%>%
    arrange(trap_ID,date_time)%>%
    group_by(trap_ID)%>%
    mutate(
      #ID gaps greater than threshold
      is_long_gap = ifelse(is.na(total_sample_minutes), FALSE, total_sample_minutes > gap_threshold_minutes),
      
      #create trapID segment where 0 is start
      trap_segment = cumsum(is_long_gap),
      
      #create decimal ID
      trap_ID_decimal = ifelse(
        trap_segment == 0,
        trap_ID,  # No long gaps yet, keep original ID
        paste0(trap_ID, ".", sprintf(as.character(trap_segment)))  #add segment
      )
    )
  
  #at this point its safe to set total_sample_minutes to 0
  #for visits that were not fishing
  data<-data%>%
    mutate(total_sample_minutes=ifelse(
      trap_status=="Not fishing",0,total_sample_minutes
    ))
  
  return(data)
}
