###############################
#function to summarize a date index
###############################
summarize_index<-function(date_index,sum.by){

  if (sum.by == "week") {
    #use ISO weeks
    index <- format(date_index, "%G-W%V")
  } else if (sum.by == "month") {
    index <- list(s.by = format(date_index, "%Y-%m"))
  } else if (sum.by == "year") {
    # For year: use mean year of the date range
    year_of_mean <- format(mean(date_index, na.rm = TRUE), "%Y")
    index <- list(s.by = rep(year_of_mean, length(date_index)))
  } else {  # day
    index <- list(s.by = format(date_index, "%Y-%m-%d"))
  }
  
  return(index)
}

