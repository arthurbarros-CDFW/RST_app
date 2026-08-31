#' @title compare_efficiency_models
#'
#' @description Models efficiency for traps using provided mark-recapture efficiency trials. 
#' 
#' @param release_data Data set with records of fish marked and 
#' released for efficiency trials.
#' 
#' @param recapture_data Data set with records of marked fish recaptured for efficiency trials.
#' 
#' @param visit_data: Input trap visit data frame.
#' 
#' @param impute_all: User input of dictating if all efficiency estimates should be imputed. 
#' Defaults to FALSE. If TRUE, will replace estimates made with actual efficiency trial data with imputed values.
#' 
#' @param survey_start Input survey season start date.
#' 
#' @param survey_end Input survey season end date.
#' 
#' @param min_sample_size Minimum number of efficiency trials required 
#' in order to fit binomial model to impute efficiency. Defaults to 10.
#' 
#' @param use_discharge User input determining if discharge should be used as a
#' covariate for modeling efficiency.
#' 
#' @param max.df.spline Maximum number of splines to fit to binomial GLM. 
#'  
#' @return 
#' 
compare_efficiency_models<-function(release_data,
                                    recapture_data,
                                    visit_data,
                                    impute_all=F,
                                    survey_start,
                                    survey_end,
                                    min_sample_size=10,
                                    use_discharge=FALSE,
                                    max.df.spline=4){
  
  ######################################
  #replicate est_efficiency to prep data for modeling
  ######################################
  
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
  
  ############################################
  #sometimes we have releases on same day but with different markers
  #so here we combine
  ############################################
  
  #find same time releases w different markers
  duplicate_releases<-trial_data%>%
    group_by(release_datetime) %>%
    summarise(
      n_release_ids = n_distinct(release_ID),
      release_ids = list(unique(release_ID))
    ) %>%
    filter(n_release_ids > 1)
  
  duplicate_releases <- trial_data %>%
    filter(release_datetime %in% duplicate_releases$release_datetime)
  
  duplicate_IDs<-unique(duplicate_releases$release_ID)
  
  duplicate_visits<-duplicate_releases%>%
    ungroup()%>%
    select(visit_datetime,trap_ID,trap_ID_decimal,trap_visit_ID,n)%>%
    group_by(visit_datetime,trap_ID,trap_ID_decimal,trap_visit_ID)%>%
    summarise(n=sum(n))
  
  duplicate_releases<-duplicate_releases%>%
    ungroup()%>%
    select(trap_ID_decimal,trap_ID,n_released,
           HrsToLastVisitAfter,HrsToFirstVisitAfter,
           release_datetime)%>%
    unique()%>%
    group_by(trap_ID_decimal,trap_ID,
             HrsToLastVisitAfter,HrsToFirstVisitAfter,
             release_datetime)%>%
    summarise(n_released=sum(n_released))
  
  duplicate_releases<-duplicate_releases%>%
    group_by(trap_ID,trap_ID_decimal,release_datetime,
             HrsToFirstVisitAfter,HrsToLastVisitAfter)%>%
    summarise(n_released=sum(n_released))
  
  duplicate_releases<-duplicate_releases%>%
    left_join(duplicate_visits)
  
  trial_data<-trial_data%>%
    filter(!release_ID %in% duplicate_IDs)
  
  trial_data<-trial_data%>%
    ungroup()%>%
    select(-release_ID)%>%
    rbind(duplicate_releases)
  
  ############################################
  
  #calculate mean recapture time weighted by number of recaps
  trial_data <- trial_data %>%
    group_by(trap_ID,trap_ID_decimal,
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
    
    #fix discharge function
    fix_discharge_field <- function(data) {
      col_idx <- which(tolower(names(data)) == "discharge")
      if (length(col_idx) > 0) {
        #rename the first match to "discharge"
        names(data)[col_idx[1]] <- "discharge"
      } else {
        warning("No 'discharge' field found in the dataset")
      }
      return(data)
    }
  
    visits <- fix_discharge_field(visits)
    
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
  
  ######################################
  #compare eff models
  ######################################
  d<-eff
  
  traps <- sort(unique(as.character(d$trap_ID_decimal)))
  ans<-list()
  ans$comparison.df <- list()
  ans$candidate_models<-list()
  ans$eff.type<-list()
  
  #set fits vector for storing outputs
  fits<-vector("list", length(traps))
  fits <- all.X <- all.ind.inside <- all.dts <- obs.data <- eff.type <- vector("list", length(traps))
  names(fits) <- traps
  names(all.X) <- traps
  names(all.dts) <- traps
  names(all.ind.inside) <- traps
  names(obs.data) <- traps
  names(eff.type) <- traps
  
  for(trap in traps){
    df<-d%>%
      filter(trap_ID_decimal==trap)
    
    #initialize model comparison
    comparison_df <- data.frame(
      Model = "No models fit",
      AIC = NA_real_,
      AICc = NA_real_,
      AICc_diff = NA_real_
    )
    
    #set "ind" which is an index of records T/F that had trials
    ind<-!is.na(df$efficiency)
    
    #set "ind.inside" which says whether a date is within trial surveys
    strt.dt<-suppressWarnings(min(df$batch_date[ind],na.rm=T))
    end.dt<-suppressWarnings(max(df$batch_date[ind],na.rm=T))
    ind.inside<-(strt.dt <= df$batch_date) & 
      (df$batch_date <= end.dt)
    inside.dates <- c(strt.dt, end.dt)
    all.ind.inside[[trap]] <- inside.dates
    
    
    #data frame for fitting (only include days with efficiency trials within season)
    tmp.df <- df[ind & ind.inside, ]
    
    #number of trials within season
    m.i <- sum(ind & ind.inside)
    
    #we are running into an issue with trials with 0 recap
    #to fix, create and effective m.i which subtracts trials with 0 recap
    zero_recaps<-nrow(df%>%filter(n_recaps==0))
    m.i=m.i-zero_recaps
    
    if( m.i == 0 ){
      ###################################################
      #If no trials, just set everything to NA
      ###################################################
      cat( paste("NO EFFICIENCY TRIALS FOR TRAP", trap, ".\n") )
      cat( paste("Catches at this trap will not be included in production estimates.\n"))
      fits[[trap]] <- NA
      all.X[[trap]] <- NA
      df$efficiency <- NA  
      obs.data[[trap]] <- NA
      eff.type[[trap]] <-0
      comparison_df <- data.frame(
        Model = "No efficiency trials",
        AIC = NA,
        AICc = NA,
        AICc_diff = NA
      )
      candidate_models <- NA
      
    }else if(m.i<min_sample_size | sum(na.omit(df$n_recaps))==0){
      ###################################################################
      #Else if too few trials (<10) for gam, use ROM+1
      ###################################################################
      cat("Fewer than 10 trials found or no recaptures.  'ROM+1' estimator used.\n")
      
      #if too few trials we assume constant efficiency over survey period
      #use Ratio of Means + 1 estimator for each trap_ID
      obs.mean <- (sum(tmp.df$n_recaps, na.rm = TRUE) + 1) / 
        (sum(tmp.df$n_released, na.rm = TRUE) + 1)
      cat(paste(trap,"'ROM+1' efficiency= ", obs.mean, "\n"))
      
      #CAMPR defaults to using the ROM for all days, missing or not,
      #so we will do the same here
      df$efficiency_imputed <- obs.mean
      
      if(impute_all==F){
        df$efficiency[!ind] <- obs.mean
      } else {
        df$efficiency<-obs.mean
      }
      
      #next we can fit a null model so we have dispersion stats, beta,
      #and covar matrix for use in later bootstrapping
      fit<-glm(efficiency~1,
               family=binomial,
               data=tmp.df,
               weights=n_released)
      
      obs.data[[trap]] <- tmp.df
      eff.type[[trap]] <-1 
      
      fits[[trap]] <- fit
      
      comparison_df <- data.frame(
        Model = "ROM+1 (constant efficiency)",
        AIC = NA,
        AICc = NA,
        AICc_diff = NA,
        type=1
      )
      candidate_models<-list(fit)
      names(candidate_models)<-comparison_df$Model
      candidate_types<-1
      names(candidate_types)<-comparison_df$Model
      
    } else {
      ####################################################################
      #there are enough efficiency trials for B-spline model temporal
      ####################################################################
      cat(paste("\n\n++++++Spline model fitting (no covariates) for trap:", trap))
      
      fits[[trap]] <- NULL
      
      #initialize discharge tracking variables
      used_discharge <- FALSE
      has_discharge <- "discharge_scaled" %in% names(tmp.df) && 
        !all(is.na(tmp.df$discharge_scaled)) && use_discharge == TRUE
      
      #1st: fit a null model
      fit <- glm( efficiency ~ 1, 
                  family=binomial, 
                  data=tmp.df, 
                  weights=n_released ) 
      fit.AIC <- AIC(fit)
      fit.AICc<-compute_aicc(fit)
      fit_type <- "null model"
      
      cat(paste("df= ", 1, ", conv= ", fit$converged, " bound= ", fit$boundary,
                " AIC= ", round(fit.AIC, 4),
                " AICc=",round(fit.AICc,4), "\n"))
      
      #added chunk from CAMPR eff_model.R that catches if model fit doesn't work
      if(!fit$converged | fit$boundary){
        cat("Constant (intercept-only) logistic model for efficiency failed. Using 'ROM+1' estimator. ")
        obs.mean<-tmp.df%>%
          summarise(ROM_eff=sum(na.omit(n_recaps+1))/sum(na.omit(n_released)+1))
        cat(paste(t,"'ROM+1' efficiency= ", obs.mean$ROM_eff, "\n"))
        
        df$efficiency <- obs.mean
        
        fits[[trap]] <- fit
        fit_type <- "rom"
        eff.type[[trap]] <- 3
        
        comparison_df <- data.frame(
          Model = "ROM+1 (constant efficiency)",
          AIC = NA,
          AICc = NA,
          AICc_diff = NA,
          type=3
        )
      } else {
        ####################################################################
        #2nd: fit models and compare
        ####################################################################
        
        #set model storage
        best_fit <- fit
        best_AIC <- fit.AIC
        best_AICc<-fit.AICc
        best_type <- "null model"
        best_used_discharge <- FALSE
        
        #track all candidate models for comparison
        candidate_models <- list()
        candidate_AIC <- c()
        candidate_AICc<-c()
        candidate_names <- c()
        candidate_types<-list()
        
        #check with AICc
        candidate_models[[1]] <- fit
        candidate_AIC[1]<-AIC(fit)
        candidate_AICc[1] <- compute_aicc(fit)  #use AICc for small sample sizes
        candidate_names[1] <- "null model"
        candidate_types[1] <- 2
        
        ####################################################################
        #fit discharge-only if applicable
        ####################################################################
        
        if(has_discharge){
          #remove rows missing discharge
          tmp.df.complete <- tmp.df[!is.na(tmp.df$discharge_scaled), ]
          
          #fit if we have enough data
          if(nrow(tmp.df.complete) >= min_sample_size){
            fit_discharge <- tryCatch(
              glm( efficiency ~ discharge_scaled, 
                   family=binomial, 
                   data=tmp.df.complete, 
                   weights=n_released ),
              error = function(e) NULL,
              warning = function(w) NULL
            )
            
            if(!is.null(fit_discharge)){
              fit_discharge.AIC<- AIC(fit_discharge)
              fit_discharge.AICc <- compute_aicc(fit_discharge)
              cat(paste("df= discharge, conv= ", fit_discharge$converged, 
                        " bound= ", fit_discharge$boundary, 
                        " AIC= ", round(fit_discharge.AIC, 4), 
                        " AICc= ", round(fit_discharge.AICc, 4), "\n"))
              
              if(fit_discharge$converged && !fit_discharge$boundary){
                candidate_models[[length(candidate_models)+1]] <- fit_discharge
                candidate_AIC <- c(candidate_AIC, fit_discharge.AIC)
                candidate_AICc <- c(candidate_AICc, fit_discharge.AICc)
                candidate_names <- c(candidate_names, "discharge_only")
                candidate_types<-c(candidate_types,3)
              }
            }
          }
        }
        
        ####################################################################
        #3rd: fit spline models with and without discharge
        ####################################################################
        #skip the quadratic df =2
        #df = 3 = cubic model (no internal knots)
        #df = 4 = cubic spline w/ 1 internal knot at median
        #df = 5 = cubic spline w/ 2 internal knots at 0.33 and 0.66 of range
        #etc. (subtract 3 from df to get number of internal knots)
        cur.df <- 3
        final_df <- 3 
        repeat{
          
          #create spline basis for all days in season
          cur.bspl <- splines::bs( df$batch_date[ind.inside],#set current b-splines
                                   df=cur.df )
          
          #subset to only days with efficiency trials
          tmp.bs <- cur.bspl[!is.na(df$efficiency[ind.inside]),] 
          
          #############################################################
          #fit temporal spline model without discharge
          #############################################################
          if(nrow(tmp.df) >= min_sample_size){
            cur.fit_temporal <- tryCatch(
              glm( efficiency ~ tmp.bs, 
                   family=binomial,
                   data=tmp.df, 
                   weights=tmp.df$n_released,
                   na.action = na.exclude, 
                   singular.ok = FALSE),
              error = function(e) NULL,
              warning = function(w) NULL
            )
            
            if(!is.null(cur.fit_temporal)){
              cur.AIC_temporal<-AIC(cur.fit_temporal)
              cur.AICc_temporal <- compute_aicc(cur.fit_temporal)
              cat(paste("Temporal df= ", cur.df, ", conv= ", cur.fit_temporal$converged, 
                        " bound= ", cur.fit_temporal$boundary, 
                        " AICc= ", round(cur.AIC_temporal, 4), "\n"))
              
              if(cur.fit_temporal$converged && !cur.fit_temporal$boundary){
                candidate_models[[length(candidate_models)+1]] <- cur.fit_temporal
                candidate_AIC <- c(candidate_AIC, cur.AIC_temporal)
                candidate_AICc <- c(candidate_AICc, cur.AICc_temporal)
                candidate_names <- c(candidate_names, paste0("temporal_df", cur.df))
                candidate_types<-c(candidate_types,4)
              }
            }
          }
          
          #############################################################
          #fit temporal spline + discharge model if available
          #############################################################
          
          if(has_discharge){
            #model with discharge and temporal spline
            tmp.df.complete <- tmp.df[!is.na(tmp.df$discharge_scaled), ]
            tmp.bs.complete <- tmp.bs[!is.na(tmp.df$discharge_scaled), ]
            
            if(nrow(tmp.df.complete)>=min_sample_size){
              cur.fit_discharge_temporal <- tryCatch(
                glm( efficiency ~ discharge_scaled + tmp.bs.complete, 
                     family=binomial,
                     data=tmp.df.complete, 
                     weights=tmp.df.complete$n_released,
                     na.action = na.exclude, 
                     singular.ok = FALSE),
                error = function(e) NULL,
                warning = function(w) NULL
              )
            }
            
            if(!is.null(cur.fit_discharge_temporal)){
              cur.AIC_discharge_temporal <- AIC(cur.fit_discharge_temporal)
              cur.AICc_discharge_temporal <- compute_aicc(cur.fit_discharge_temporal)
              cat(paste("Discharge+Spline df= ", cur.df, ", conv= ", cur.fit_discharge_temporal$converged, 
                        " bound= ", cur.fit_discharge_temporal$boundary, 
                        " AIC= ", round(cur.AIC_discharge_temporal, 4),
                        " AICc= ", round(cur.AICc_discharge_temporal, 4), "\n"))
              
              if(cur.fit_discharge_temporal$converged && !cur.fit_discharge_temporal$boundary){
                candidate_models[[length(candidate_models)+1]] <- cur.fit_discharge_temporal
                candidate_AIC <- c(candidate_AIC, cur.AIC_discharge_temporal)
                candidate_AICc <- c(candidate_AICc, cur.AICc_discharge_temporal)
                candidate_names <- c(candidate_names, paste0("discharge_temporal_df", cur.df))
                candidate_types<-c(candidate_types,5)
              }
            }
          }
          
          #check repeat stop conditions
          if(cur.df > max.df.spline){
            break
          }
          cur.df <- cur.df + 1
        }
        
        ###############################################################
        #select best model based on AIC
        ###############################################################
        
        if(length(candidate_models)>1){
          #find model with min AIC
          best_idx<-which.min(candidate_AICc)
          best_fit<-candidate_models[[best_idx]]
          best_AIC<-candidate_AIC[[best_idx]]
          best_AICc<-candidate_AICc[[best_idx]]
          best_name <- candidate_names[best_idx]
          best_type <- candidate_types[best_idx]
          
          best_used_discharge <- any(grepl("discharge", names(coef(best_fit))))
          
          cat("\n========================================\n")
          cat("Best model for trap:", trap, "\n")
          cat("Model type:", best_name, "\n")
          cat("AICc:", round(best_AIC, 4), "\n")
          
          #compare top models
          comparison_df <- data.frame(
            Model = candidate_names,
            AIC=round(candidate_AIC,4),
            AICc = round(candidate_AICc, 4),
            AICc_diff = round(candidate_AICc - min(candidate_AICc), 2)
          )
          comparison_df <- comparison_df[order(comparison_df$AICc), ]
        } else {
          #just in case only one candidate
          best_fit <- candidate_models[[1]]
          best_AIC <- candidate_AIC[1]
          best_AICc <- candidate_AICc[1]
          best_type <- candidate_names[1]
          best_used_discharge <- FALSE
          comparison_df <- data.frame(
            Model = candidate_names,
            AIC = round(candidate_AIC, 4),
            AICc = round(candidate_AICc, 4),
            AICc_diff = 0
          )
        }
        
        if(exists("best_fit") && !is.null(best_fit)) {
          fits[[trap]] <- best_fit
        } else {
          # Fallback to the null model
          fits[[trap]] <- fit
        }
        
        #save info about model selection
        if(exists("candidate_models") && length(candidate_models) > 0){
          attr(fits[[trap]], "candidate_
               ") <- data.frame(
                 Model = candidate_names,
                 AIC = round(candidate_AIC, 4),
                 AICc = round(candidate_AICc, 4)
               )
          attr(fits[[trap]], "model_selected") <- best_type
          attr(fits[[trap]], "used_discharge") <- used_discharge
        }
        
        names(candidate_models)<-candidate_names
        names(candidate_types)<-candidate_names
      }
      #save the raw efficiency data.  
      obs.data[[trap]] <- tmp.df
      eff.type[[trap]] <- 4
    }
    
    if(is.na(fits[trap])){  
      ans$comparison.df[[trap]]<-NA
      
    } else {
      
      ans$comparison.df[[trap]]<-comparison_df
      ans$candidate_models[[trap]]<-candidate_models
      ans$eff.type[[trap]]<-candidate_types
    }
    
  }
  
  ans$eff<-eff
  return(ans)

}
