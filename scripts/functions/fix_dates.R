#' @title fix_dates
#'
#' @description Converts date fields of input data frame to POSIXct with specified time zone and format (%Y-%m-%dT%H:%M:%SZ). Called for by the est_catch.R and est_efficiency functions.
#' 
#' @param df Input data frame with date field needed to be cleaned/formatted.
#' 
#' @param time_field Field with visit times to be formatted, needs to be titled "visit_time".
#' 
#' @return df Outputs data frame with formatted time field.
#' 
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
