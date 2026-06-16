#' @title summarize_passage
#'
#' @description Aggregate passage estimates by sum.by time period. Takes the following steps:
#' 1) Averages trap passage estimates over batch_date.
#' 2) Creates date index using summarize_index().
#' 3) Summarize passage by sum.by time period.
#' 
#' @param passage_data  Date value to be summarized to grouping unit.
#' 
#' @param sum.by Time grouping variable to summarize passage data by (default: "week").
#' 
#' @return n Data frame with sum.by value, passage estimates, dates, and estimate of the percent of imputed catch for the grouping unit.
#' 
summarize_passage<-function(passage_data,sum.by){
  
  index <- list(batch_date=format( passage_data$batch_date, "%Y-%m-%d" ))
  
  #average over batch_date across traps
  n <- c(tapply( passage_data$passage, index, FUN=mean, na.rm=T ))
  dt <- c(tapply( passage_data$batch_date, index, #get minimum date per day
                  FUN=function(x, narm){ min(x, na.rm=narm) }, narm=T ))
  class(dt) <- class(passage_data$batch_date)  # the tapply above strips class information, so here make this back into a POSIX date
  attr(dt, "tzone") <- attr(passage_data$batch_date, "tzone")  
  
  attr(dt, "jday") <- attr(passage_data, "jday") 
  p.imp <- c(tapply( passage_data$p.c.imputed, index,
                     FUN=mean, na.rm=T ))   # p.c.imputed is the var with % per day
  
  #now summarize passage by sum.by time period
  index<-summarize_index(date_index = dt,sum.by)
  n <- c(tapply( n, index, FUN=sum, na.rm=T )) #sum n by sum.by time groups
  dt <- c(tapply( dt, index, #get minimum date by group
                  FUN=function(x, narm){ min(x, na.rm=narm) }, narm=T )) 
  p.imp <- c(tapply( p.imp, index, FUN=mean, na.rm=T )) #get average %imputed per groups
  class(dt) <- class(passage_data$batch_date)
  
  #return
  n <- data.frame( s.by=names(n),
                   passage=n,
                   date=dt,
                   pct.imputed.catch=p.imp,
                   stringsAsFactors=F, row.names=NULL)
  
  n
}
