###############################
#function sets dates to target time zone and format
###############################
fix_dates<-function(df,time_field){
  if(!time_field %in% names(df)){
    stop("The specified time_field '", time_field, "' does not exist in the dataframe")
  }
  time_values <- df[[time_field]]
  df<-df%>%
    mutate(date_time=as.POSIXct(time_values,
                               format = "%Y-%m-%dT%H:%M:%SZ",
                               tz = time_zone))
  
  return(df)   
}
