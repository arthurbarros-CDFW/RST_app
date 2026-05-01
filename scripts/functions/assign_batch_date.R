###############################
#function assigns "batch" dates to date data
#used to link sampling visits and catch data to efficiency trials
###############################
assign_batch_date<-function(df,time_field){
  if(!time_field %in% names(df)){
    stop("The specified time_field '", time_field, "' does not exist in the dataframe")
  }
  
  time_zone<-"America/Los_Angeles"
  #batch date using a 4am rule
  cut_time <- "04:00:00"
  
  #CAMPR makes this little "safe.ifelse" function because ifelse drops the posixct format for some reason
  safe.ifelse <- function(cond, yes, no) {
    structure(ifelse(cond, yes, no), class = class(yes)) 
  }
  
  time_values <- df[[time_field]]
  
  if (!inherits(time_values, "POSIXct")) {
    stop("The specified time_field '", time_field, "' must be of class POSIXct")
  }
  df$batch_date <- as.Date(
    safe.ifelse(
      format(time_values, "%H:%M:%S") <= cut_time,
      time_values - 24 * 60 * 60,  #subtract 1 day if before cutoff
      time_values                   #keep same day if after cutoff
    ),
    tz = time_zone
  )
  
  return(df)
}