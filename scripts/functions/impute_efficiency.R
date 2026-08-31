#' @title impute_efficiency
#'
#' @description 
#' 
#' @param 
#' 
#' @param 
#' 
#' @param 
#' 
#' @param 
#' 
#' @return 
#' 

impute_efficiency<-function(efficiency_data,
                           selected_models,
                           model_names,
                           impute_all=F){
  
  #set variables
  d<-efficiency_data
  models<-selected_models
  traps <- sort(unique(as.character(d$trap_ID_decimal)))
  ans<-list()
  ans$results<-NULL
  ans$models<-list()
  warnings_collected <- list()
  
  all.X <- all.ind.inside <- all.dts <- obs.data <- eff.type <- vector("list", length(traps))
  names(all.X) <- traps
  names(all.dts) <- traps
  names(all.ind.inside) <- traps
  names(obs.data) <- traps
  names(eff.type) <- traps
  
  #we are going to do efficiency for each trap_ID separately
  for(trap in traps){
    df<-d%>%
      filter(trap_ID_decimal==trap)
    
    fit<-models[[trap]]
    fit_name<-model_names[[trap]]
    
    used_discharge <- any(grepl("discharge", names(coef(fit))))
    
    #set "ind" which is an index of records T/F that had trials
    ind<-!is.na(df$efficiency)
    
    #set "ind.inside" which says whether a date is within trial surveys
    strt.dt<-suppressWarnings(min(df$batch_date[ind],na.rm=T))
    end.dt<-suppressWarnings(max(df$batch_date[ind],na.rm=T))
    ind.inside<-(strt.dt <= df$batch_date) & 
      (df$batch_date <= end.dt)
    inside.dates <- c(strt.dt, end.dt)
    all.ind.inside[[trap]] <- inside.dates
    
    tmp.df <- df[ind & ind.inside, ]
    
    n_coef <- length(coef(fit))
    
    if( n_coef <= 1 ){
      #intercept only model
      X <- matrix( 1, sum(ind.inside), 1)
    } else {
      
      bspl_final <- splines::bs( df$batch_date[ind.inside], df = ifelse(n_coef > 2, n_coef - 1, 1) )
      
      if(used_discharge){
        #discharge + spline model
        discharge_full <- df$discharge_scaled[ind.inside]
        discharge_full[is.na(discharge_full)] <- 0
        X <- cbind(1, discharge_full, bspl_final)
      } else if(grepl("discharge_only", fit_name)){
        #discharge-only model
        discharge_full <- df$discharge_scaled[ind.inside]
        discharge_full[is.na(discharge_full)] <- 0
        X <- cbind(1, discharge_full)
      } else {
        #temporal spline model
        X <- cbind(1, bspl_final)
      }
      
      #check X has correct number of columns matching coefficients
      if(ncol(X) != n_coef){
        #if mismatch, try to rebuild
        if(used_discharge){
          #rebuild with proper dimensions
          n_spline_cols <- n_coef - 2  # intercept + discharge
          if(n_spline_cols > 0){
            bspl_final <- splines::bs( df$batch_date[ind.inside], df = n_spline_cols )
            discharge_full <- df$discharge_scaled[ind.inside]
            discharge_full[is.na(discharge_full)] <- 0 #sets missing discharge to avg (scaled)
            X <- cbind(1, discharge_full, bspl_final)
            
            colnames(X)<-names(coef(fit))
          } else {
            X <- cbind(1, discharge_full)
          }
        } else {
          #temporal model
          n_spline_cols <- n_coef - 1
          bspl_final <- splines::bs( df$batch_date[ind.inside], df = n_spline_cols )
          X <- cbind(1, bspl_final)
        }
      }
    }
    
    #make predictions
    pred <- X %*% coef(fit)
    pred <- 1 / (1 + exp(-pred))
    
    #store observed discharge range for extrapolation checking
    if(used_discharge){
      observed_discharge_range <- range(tmp.df$discharge_scaled, na.rm=TRUE)
      attr(fit, "observed_discharge_range") <- observed_discharge_range
      cat("Observed discharge range for trap", trap, ":", 
          round(observed_discharge_range[1], 2), "to", 
          round(observed_discharge_range[2], 2), "(scaled)\n")
    }
    
    #check for discharge extrapolation if discharge was used
    if(used_discharge){
      observed_range <- attr(fit, "observed_discharge_range")
      if(!is.null(observed_range)){
        #get discharge values being predicted for
        pred_discharge <- df$discharge_scaled[ind.inside]
        
        #find which days are outside observed range
        outside_low <- pred_discharge < observed_range[1] & !is.na(pred_discharge)
        outside_high <- pred_discharge > observed_range[2] & !is.na(pred_discharge)
        outside_range <- outside_low | outside_high
        
        if(any(outside_range, na.rm=TRUE)){
          #get dates for outside-range predictions
          outside_dates <- df$batch_date[ind.inside][outside_range]
          outside_dates <- outside_dates[!is.na(outside_dates)]
          
          warning_msg <- paste0(
            "EXTRAPOLATION WARNING for trap: ", trap, "\n",
            "  Predicting efficiency for ", sum(outside_range, na.rm=TRUE), 
            " days outside observed discharge range\n",
            "  Observed range: ", round(observed_range[1], 2), " to ", 
            round(observed_range[2], 2), " (scaled)\n",
            "  Prediction range: ", round(min(pred_discharge, na.rm=TRUE), 2), 
            " to ", round(max(pred_discharge, na.rm=TRUE), 2), " (scaled)\n",
            "  Percent of predictions outside range: ", 
            round(sum(outside_range, na.rm=TRUE) / sum(!is.na(pred_discharge)) * 100, 1), "%"
          )
          
          # Store warning
          warnings_collected[[length(warnings_collected) + 1]] <- warning_msg
          
          # Still print to console for debugging
          cat("\n⚠️ ", warning_msg, "\n")
          cat("========================================\n")
        }
      }
    }
    
    #replace all days efficiency (within trial days) with predicted
    df$efficiency_imputed[ind.inside] <- pred
    
    #use the mean of predictions for all dates outside trial period
    mean.p<-mean(pred,na.rm=T)
    df$efficiency_imputed[!ind.inside]<-mean.p
    
    #impute all or not
    if(impute_all==F & !is.null(fit)){
      df$effiency_imputed <- factor( !ind, levels=c(T,F), labels=c(TRUE, FALSE))
      df$efficiency<-ifelse(df$effiency_imputed==TRUE,
                            df$efficiency_imputed,df$efficiency)
      #if original trial efficiency is 0, impute, this is not done in
      #campR because everything is imputed
      #we have to do this because eff of 0 makes inf catch
      df <- df %>%
        mutate(
          effiency_imputed = ifelse(efficiency == 0, TRUE, 
                                    as.character(effiency_imputed)),
          efficiency = ifelse(efficiency == 0, efficiency_imputed, efficiency)
        )
    }else{
      df$effiency_imputed<-TRUE
      df$efficiency<-df$efficiency_imputed
    }
    
    if(!is.null(fit)){
      #base columns that always exist
      base_cols <- c("trap_ID_decimal", "batch_date", "n_released", "n_recaps", 
                     "efficiency_imputed","unimputed_efficiency", "efficiency")
      
      #add discharge columns if exist
      if("discharge" %in% names(df) && "discharge_scaled" %in% names(df)) {
        select_cols <- c("trap_ID_decimal", "batch_date", "n_released", "n_recaps", 
                         "discharge", "discharge_scaled", "unimputed_efficiency",
                         "efficiency_imputed", "efficiency")
      } else {
        select_cols <- base_cols
      }
      
      df <- df[, select_cols]
    }
    
    
    if(!is.null(fit)){
      df <- df[, select_cols]
    }
    
    if(is.null(fit)){  
      #save the fit for bootstrapping.
      ans$models[[trap]] <- NA
      ans$all.X[[trap]] <- NA
      ans$all.ind.inside <- all.ind.inside
      ans$all.dts[[trap]] <- df$batch_date[ind.inside]
      ans$obs.data[[trap]] <- tmp.df
      ans$eff.type <- eff.type
      ans$comparison.df[[trap]]<-NA
      ans$warnings[[trap]]<-warnings_collected
      
      #base columns that always exist
      base_cols <- c("trap_ID", "trap_ID_decimal", "batch_date", "n_released", "n_recaps", "efficiency_imputed")
      
      #add discharge columns if they exist
      if("discharge" %in% names(df) && "discharge_scaled" %in% names(df)) {
        select_cols <- c("trap_ID", "trap_ID_decimal", "batch_date", "n_released", "n_recaps", 
                         "discharge", "discharge_scaled", "effiency_imputed")
      } else {
        select_cols <- base_cols
      }
      
      df <- df[, select_cols]
      df <- df %>% mutate(effiency_imputed = NA,
                          efficiency = NA,
                          enhanced.eff = NA)
      ans$results <- ans$results %>% rbind(df)
      
    } else {
      #save the fit for bootstrapping.
      ans$models[[trap]] <- fit
      ans$all.X[[trap]] <- X
      ans$all.ind.inside <- all.ind.inside
      ans$all.dts[[trap]] <- df$batch_date[ind.inside]
      ans$obs.data[[trap]] <- tmp.df
      ans$warnings[[trap]]<-warnings_collected
      
      ans$results <- ans$results %>% rbind(df)
    }

  }
  
  return(ans)
}
