#' @title summarize_index
#'
#' @description Convert dates to grouping units (day, week, month, year). Called by summarize_passage().
#' 
#' @param date_index  Date value to be summarized to grouping unit.
#' 
#' @param sum.by Time grouping variable to summarize passage data by (default: "week").
#' 
#' @return index Assigned date grouping value to date records.
#' 
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

