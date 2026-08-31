#' @title est_passage
#'
#' @description Coordinate all modules to produce final passage estimates. Takes the following steps:
#' 1) Run est_catch() to calculate catch for available data and impute catch for missing periods.
#' 2) Run est_efficiency() to calculate efficiency for available trials and impute for missing periods.
#' 3) Estimate passage and related uncertainty.
#' 
#' @param catch Records of fish caught during each trap visit.
#' 
#' @param visits  Records of all trap operation and sampling events.
#' 
#' @param release Records of fish marked and released for efficiency trials.
#' 
#' @param recapture Records of marked fish recaptured for efficiency trials.
#' 
#' @param summarize.by Time grouping variable to summarize passage data by (default: "week").
#' 
#' @param impute_all Determine whether or not to impute all efficiency values, even those observed via efficiency trials (default: FALSE).
#' 
#' @param bootstrap Whether or not to conduct bootstrapping for uncertainty (default: TRUE).
#' 
#' @param survey_start Input survey season start date.
#' 
#' @param survey_end Input survey season end date.
#' 
#' @param target_species Species of interest, currently limited to Chinook Salmon.
#' 
#' @param target_run Run of interest, currently limited to Fall.
#' 
#' @param file.name File name for saved outputs, which are currently commented out (default: "test").
#' 
#' @return results List with estimates of passage, catch, and efficiency.

est_passage<-function(catch.results,
                      eff.results,
                      summarize.by="day",
                      bootstrap=T,survey_start,survey_end,
                      catch.fits,
                      eff.fits,
                      eff.types,
                      target_species,
                      target_run=NA,
                      use_discharge=F,
                      min_sample_size=10){
  
  #set sum.by ifelse to change day to jday
  #this is so user can just put day but code sees julian day for easier editing
  #this may be messy
  sum.by=ifelse(summarize.by=="day","jday",summarize.by)
  
  #pull relvant objects from catch.results
  catch.fits=catch.results$models
  catch.X.miss=catch.results$X.miss
  catch.gapLens=catch.results$gaps
  catch.bDates.miss=catch.results$batchDate.for.missings
  catch.results<-catch.results$results
  
  #pull relevant objects from eff.results
  eff.fits=eff.results$models
  eff.X=eff.results$all.X
  eff.ind.inside=eff.results$all.ind.inside
  eff.X.dates=eff.results$all.dts
  eff.X.obs.data=eff.results$obs.data
  eff.results=eff.results$results
  
  catch.results<-catch.results%>%
    select(trap_ID_decimal,batch_date,trap_status,
           total_sample_minutes,total_catch,catch_imputed)
  
  if(use_discharge){
    eff.results<-eff.results%>%
      select(trap_ID_decimal,batch_date,discharge,discharge_scaled,
             efficiency_imputed,efficiency,
             unimputed_efficiency)
  }else{
    eff.results<-eff.results%>%
      select(trap_ID_decimal,batch_date,
             efficiency_imputed,efficiency,
             unimputed_efficiency)
  }

  
  #filter catch.results by start and end dates, first time we do this
  catch.results<-catch.results%>%
    filter(batch_date<=survey_end & batch_date>=survey_start)
  
  pass_data<-catch.results%>%
    left_join(eff.results)
  
  pass_data$passage<-NA
  pass_data<-pass_data%>%
    mutate(passage=total_catch/efficiency)
  pass_data$passage<-round(pass_data$passage)
  
  #estimate percent imputed catch per day (for most cases this will be
  #for just one sampling event as we dont have more than one visit per day)
  pass_data <- pass_data %>%
    group_by(batch_date, trap_ID_decimal) %>%
    mutate(p.c.imputed = mean(catch_imputed),
           p.e.imputed = mean(efficiency_imputed)) %>%
    ungroup()

  #note: may be able to remove below up until mark
  ###############
  #summarize passage data replicating summarize_passage.r

  #we're going to do things differently. CAMPR summarizes depending on what the user inputs.
  #I want to just make all the summaries, and than report based on what user wants
  pass_data$jday<-as.POSIXlt(pass_data$batch_date)$yday + 1 #julian date
  pass_data$week<-as.numeric(strftime(pass_data$batch_date, format = "%W"))
  pass_data$month<-as.numeric(strftime(pass_data$batch_date,format = "%m"))
  pass_data$year<-as.numeric(strftime(pass_data$batch_date,format = "%Y"))
  
  total_hours_sampled<-sum(pass_data$total_sample_minutes)/60
  
  fit_traps<-unique(pass_data$trap_ID_decimal)

  #plot passage estimates
  p_passage<-ggplot()+
    geom_point(data=pass_data,aes(x=batch_date,
                                  y=passage,
                                  color=catch_imputed))+
    theme_bw()+
    facet_wrap(.~trap_ID_decimal)
  
  #next we should be passing to bootstrapping passage_boot.R
  n<-passage_boot(passage_data=pass_data,
                 sum.by,catch.fits,
                  catch.X.miss=catch.X.miss,
                  catch.gapLens=catch.gapLens,
                  catch.bDates.miss=catch.bDates.miss,
                  eff.fits=eff.fits,
                  eff.X=eff.X,
                  eff.ind.inside=eff.ind.inside,
                  eff.X.dates=eff.X.dates,
                  eff.X.obs.data=eff.X.obs.data,
                 eff.types=eff.types,
                  survey_start,survey_end,
                  R=100,conf=0.95,ci=T,
                 use_discharge=use_discharge,
                 min_sample_size = min_sample_size)
  
  
  #plot breaking showing passage estimate and CI broken down by summarize.by val
  p_passage_boot<-ggplot(data=n,aes(x=s.by,y=passage))+
    geom_point()+
    geom_errorbar(aes(ymin = lower.95, 
                      ymax = upper.95), 
                  width = 0.2)+
    labs(x=summarize.by)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  p_passage_boot
  
  #ggsave(paste("outputs/boot_",Sys.Date(),"_",file.name,".png",sep=""),
   #      p_passage_boot,scale=2)
  
  #write.csv(n,paste("outputs/passage_output_",Sys.Date(),"_",file.name,".csv",sep=""),
    #        row.names = F)
  
  results<-list("p_passage"=p_passage_boot,
                "passage_output"=n)
  return(results)
}
