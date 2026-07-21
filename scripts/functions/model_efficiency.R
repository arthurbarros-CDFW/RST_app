#' @title model_efficiency
#'
#' @description Fits Binomial GLM to efficiency trial data and imputes efficiency for missing periods. Called by est_efficiency().
#' Takes the following steps:
#' 1) Fits a null model to efficiency data for each trap sampling segment.
#' 2) Calculates AIC for model fitting.
#' 3) If enough real efficiency trial data available, increase GLM spline degrees of freedom until model fit is good enough
#'  (ie model doesn’t converge or AIC increases etc.). If too few trials (<10) utilizes ration of means +1 assuming constant efficiency.
#' 4) Once model is fit, impute efficiency by looping over batch date sampling gaps.
#' 
#' @param efficiency_data Input of efficiency data frame with periods of missing efficiency data.
#' 
#' @param max.df.spline Maximum number of splines to fit to binomial GLM.
#' 
#' @param impute_all User input of dictating if all efficiency estimates should be imputed. Defaults to FALSE. If TRUE, will replace estimates made with actual
#' efficiency trial data with imputed values.
#' 
#' @param min_sample_size  Minimum number of efficiency trials required in order to fit binomial model to impute efficiency. Defaults to 10.
#' 
#' @return ans List containing efficiency model fit results and efficiency data with imputed efficiency values.
#' Also contains X.miss, model matrix with intercept and predictor terms for generating efficiency predictions.

min_sample_size=10 #default minimum sample size for efficiency trials in CAMPR
#if we have 10 efficiency trials or more we will model efficiency for missing dates
#otherwise use a constant efficiency
#NOTE: The variable names and structure for a lot of this is written based on Trent's code, I'd like to work through it and make it more legible once done

model_efficiency<-function(efficiency_data,
                           max.df.spline=4,
                           impute_all=F,
                           min_sample_size=10,
                           use_discharge=FALSE,
                           force_discharge=FALSE){
  
  #set variables
  d<-efficiency_data
  traps <- sort(unique(as.character(d$trap_ID_decimal)))
  ans<-list()
  ans$results<-NULL
  ans$models<-list()
  
  #set fits vector for storing outputs
  fits<-vector("list", length(traps))
  fits <- all.X <- all.ind.inside <- all.dts <- obs.data <- eff.type <- vector("list", length(traps))
  names(fits) <- traps
  names(all.X) <- traps
  names(all.dts) <- traps
  names(all.ind.inside) <- traps
  names(obs.data) <- traps
  names(eff.type) <- traps
  
  #we are going to do efficiency for each trap_ID separately
  for(trap in traps){
    df<-d%>%
      filter(trap_ID_decimal==trap)
    
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
      eff.type[[trap]] <- 1
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
      eff.type[[trap]] <- 2
      
      fits[[trap]] <- fit
      
      if( length(coef(fits[[trap]])) == 1 ){
        #pred <- matrix( coef(fit), sum(ind.inside), 1 )
        X <- matrix( 1, sum(ind.inside), 1)
      }
      
      #   ---- Save X, and the dates at which we predict, for bootstrapping.
      all.X[[trap]] <- X   
      all.dts[[trap]] <- df$batch_date[ind.inside] 
    } else {
      ####################################################################
      #there are enough efficiency trials for B-spline model temporal
      ####################################################################
      cat(paste("\n\n++++++Spline model fitting (no covariates) for trap:", trap))
      
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
      fit_type <- "null"
      
      cat(paste("df= ", 1, ", conv= ", fit$converged, " bound= ", fit$boundary, " AIC= ", round(fit.AIC, 4), "\n"))
      
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
      } else {
        ####################################################################
        #2nd: fit models and compare
        ####################################################################
        
        #set model storage
        best_fit <- fit
        best_AIC <- fit.AIC
        best_type <- "null"
        best_used_discharge <- FALSE
        
        #track all candidate models for comparison
        candidate_models <- list()
        candidate_AICs <- c()
        candidate_names <- c()
        
        #check with AICc
        candidate_models[[1]] <- fit
        candidate_AICs[1] <- compute_aicc(fit)  # Use AICc for small sample sizes
        candidate_names[1] <- "null"
        
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
              fit_discharge.AIC <- compute_aicc(fit_discharge)
              cat(paste("df= discharge, conv= ", fit_discharge$converged, 
                        " bound= ", fit_discharge$boundary, 
                        " AICc= ", round(fit_discharge.AIC, 4), "\n"))
              
              if(fit_discharge$converged && !fit_discharge$boundary){
                candidate_models[[length(candidate_models)+1]] <- fit_discharge
                candidate_AICs <- c(candidate_AICs, fit_discharge.AIC)
                candidate_names <- c(candidate_names, "discharge_only")
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
              cur.AIC_temporal <- compute_aicc(cur.fit_temporal)
              cat(paste("Temporal df= ", cur.df, ", conv= ", cur.fit_temporal$converged, 
                        " bound= ", cur.fit_temporal$boundary, 
                        " AICc= ", round(cur.AIC_temporal, 4), "\n"))
              
              if(cur.fit_temporal$converged && !cur.fit_temporal$boundary){
                candidate_models[[length(candidate_models)+1]] <- cur.fit_temporal
                candidate_AICs <- c(candidate_AICs, cur.AIC_temporal)
                candidate_names <- c(candidate_names, paste0("temporal_df", cur.df))
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
              cur.fit_discharge_spline <- tryCatch(
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
            
            if(!is.null(cur.fit_discharge_spline)){
              cur.AIC_discharge_spline <- compute_aicc(cur.fit_discharge_spline)
              cat(paste("Discharge+Spline df= ", cur.df, ", conv= ", cur.fit_discharge_spline$converged, 
                        " bound= ", cur.fit_discharge_spline$boundary, 
                        " AICc= ", round(cur.AIC_discharge_spline, 4), "\n"))
              
              if(cur.fit_discharge_spline$converged && !cur.fit_discharge_spline$boundary){
                candidate_models[[length(candidate_models)+1]] <- cur.fit_discharge_spline
                candidate_AICs <- c(candidate_AICs, cur.AIC_discharge_spline)
                candidate_names <- c(candidate_names, paste0("discharge_spline_df", cur.df))
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
          best_idx<-which.min(candidate_AICs)
          best_fit<-candidate_models[[best_idx]]
          best_AIC<-candidate_AICs[[best_idx]]
          best_type <- candidate_names[best_idx]
          
          best_used_discharge <- any(grepl("discharge", names(coef(best_fit))))
          
          cat("\n========================================\n")
          cat("Best model selected for trap:", trap, "\n")
          cat("Model type:", best_type, "\n")
          cat("AICc:", round(best_AIC, 4), "\n")
          
          #compare top models
          comparison_df <- data.frame(
            Model = candidate_names,
            AICc = round(candidate_AICs, 4),
            AICc_diff = round(candidate_AICs - min(candidate_AICs), 2)
          )
          comparison_df <- comparison_df[order(comparison_df$AICc), ]
          print(comparison_df)
          cat("========================================\n")
        } else {
          #just in case only one candidate
          best_fit <- candidate_models[[1]]
          best_AIC <- candidate_AICs[1]
          best_type <- candidate_names[1]
          best_used_discharge <- FALSE
        }
        
        #use best fit for predictions
        fit <- best_fit
        fit_type <- best_type
        used_discharge <- best_used_discharge
        
        #store observed discharge range for extrapolation checking
        if(used_discharge && has_discharge){
          observed_discharge_range <- range(tmp.df$discharge_scaled, na.rm=TRUE)
          attr(fit, "observed_discharge_range") <- observed_discharge_range
          cat("Observed discharge range for trap", trap, ":", 
              round(observed_discharge_range[1], 2), "to", 
              round(observed_discharge_range[2], 2), "(scaled)\n")
        }
        
        #################################################
        #make a design matrix for ease in calculating predictions.
        #################################################
        n_coef <- length(coef(fit))
        if( n_coef <= 1 ){
          #intercept only model
          X <- matrix( 1, sum(ind.inside), 1)
        } else {
          
          bspl_final <- splines::bs( df$batch_date[ind.inside], df = ifelse(n_coef > 2, n_coef - 1, 1) )
          
          if(used_discharge && has_discharge){
            #discharge + spline model
            discharge_full <- df$discharge_scaled[ind.inside]
            discharge_full[is.na(discharge_full)] <- 0
            X <- cbind(1, discharge_full, bspl_final)
          } else if(grepl("discharge_only", best_type)){
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
            if(used_discharge && has_discharge){
              #rebuild with proper dimensions
              n_spline_cols <- n_coef - 2  # intercept + discharge
              if(n_spline_cols > 0){
                bspl_final <- splines::bs( df$batch_date[ind.inside], df = n_spline_cols )
                discharge_full <- df$discharge_scaled[ind.inside]
                discharge_full[is.na(discharge_full)] <- 0
                X <- cbind(1, discharge_full, bspl_final)
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
        
        #check for discharge extrapolation if discharge was used
        if(used_discharge && has_discharge){
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
              
              cat("\n⚠️  EXTRAPOLATION WARNING for trap:", trap, "\n")
              cat("  Predicting efficiency for", sum(outside_range, na.rm=TRUE), 
                  "days outside observed discharge range\n")
              cat("  Observed range:", round(observed_range[1], 2), "to", 
                  round(observed_range[2], 2), "(scaled)\n")
              cat("  Prediction range:", round(min(pred_discharge, na.rm=TRUE), 2), 
                  "to", round(max(pred_discharge, na.rm=TRUE), 2), "(scaled)\n")
              
              #show first few dates with extrapolation
              if(length(outside_dates) > 0){
                cat("  First 5 extrapolated dates:", 
                    format(head(outside_dates, 5), "%Y-%m-%d"), "\n")
              }
              
              #store extrapolation info in the fit object for later reference
              attr(fit, "extrapolation_dates") <- outside_dates
              attr(fit, "extrapolation_count") <- sum(outside_range, na.rm=TRUE)
              attr(fit, "extrapolation_percent") <- sum(outside_range, na.rm=TRUE) / 
                sum(!is.na(pred_discharge)) * 100
              cat("  Percent of predictions outside range:", 
                  round(attr(fit, "extrapolation_percent"), 1), "%\n")
              cat("========================================\n")
            }
          }
        }
        
        #replace all days efficiency (within trial days) with predicted
        df$efficiency_imputed[ind.inside] <- pred
        
        #use the mean of predictions for all dates outside trial period
        mean.p<-mean(pred,na.rm=T)
        df$efficiency_imputed[!ind.inside]<-mean.p
        
        #note: I believe the following is a bug and have commented it out
        #it is to reduce the dataset after expanding dates, but the
        #expansion was only done in the enhanced efficiency model chunk
        #df <- df[df$batchDate >= min.date.p & df$batchDate <= max.date.p,]
        
        #save used fit for trap
        fits[[trap]] <- fit
        
        #save info about model selection
        if(exists("candidate_models") && length(candidate_models) > 0){
          attr(fits[[trap]], "candidate_models") <- data.frame(
            Model = candidate_names,
            AIC = round(candidate_AICs, 4)
          )
          attr(fits[[trap]], "model_selected") <- best_type
          attr(fits[[trap]], "used_discharge") <- used_discharge
        }
      }
      #save the raw efficiency data.  
      obs.data[[trap]] <- tmp.df
      eff.type[[trap]] <- 4
    }
    #this is legacy
    #the CAMPR script leaves it uncommented, assuming that 
    #we will always use imputed values for all days.
    #I will want to ensure that this is an option, user can either
    #impute all days OR impute just days missing trials
    #i don't like that we reuse the ind value here that we 
    #set earlier
    if(impute_all==F & !is.na(fits[trap])){
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
    
    if(!is.na(fits[trap])){
      #base columns that always exist
      base_cols <- c("trap_ID", "batch_date", "n_released", "n_recaps", 
                     "efficiency_imputed","unimputed_efficiency", "efficiency", "trap_ID_decimal")
      
      #add discharge columns if exist
      if("discharge" %in% names(df) && "discharge_scaled" %in% names(df)) {
        select_cols <- c("trap_ID", "batch_date", "n_released", "n_recaps", 
                         "discharge", "discharge_scaled", 
                         "efficiency_imputed", "efficiency", "trap_ID_decimal")
      } else {
        select_cols <- base_cols
      }
      
      df <- df[, select_cols]
    }
    
    df$trap_ID_decimal <- trap
    #df$effiency_imputed<-as.logical(df$effiency_imputed)
    
    if(!is.na(fits[trap])){
      df <- df[, select_cols]
    }
    
    if(is.na(fits[trap])){  
      #save the fit for bootstrapping.
      ans$models[[trap]] <- NA
      ans$all.X[[trap]] <- NA
      ans$all.ind.inside <- all.ind.inside
      ans$all.dts[[trap]] <- df$batch_date[ind.inside]
      ans$obs.data[[trap]] <- tmp.df
      ans$eff.type <- eff.type
      
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
      ans$eff.type <- eff.type
      
      ans$results <- ans$results %>% rbind(df)
    }
    
  }
  
  return(ans)
}
