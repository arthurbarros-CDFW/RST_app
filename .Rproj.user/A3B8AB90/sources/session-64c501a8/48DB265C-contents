###############################
#function to model and impute efficiency
###############################
min_sample_size=10 #default minimum sample size for efficiency trials in CAMPR
#if we have 10 efficiency trials or more we will model efficiency for missing dates
#otherwise use a constant efficiency
#NOTE: The variable names and structure for a lot of this is written based on Trent's code, I'd like to work through it and make it more legible once done

model_efficiency<-function(efficiency_data,
                           max.df.spline=4,
                           impute_all=F,
                           min_sample_size=10){
  
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
  
  ###############################
  #the first thing eff_model.R does is look to see if we already have
  #betas from prior efficiency models. If we do not, than we know it
  #is a new trap/subsite meaning no enhanced efficiency model for it,
  #and no data in the env covariate database.
  #We will figure this out later to determine if it's necessary.
  ###############################
  
  #we are going to do efficiency for each trap_ID separately
  for(trap in traps){
    df<-d%>%
      filter(trap_ID_decimal==trap)
    
    #set "ind" which is the list of dates T/F that had trials
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
      fit<-glm(n_recaps/n_released~1,
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
      
      #1st: fit a null model
      fit <- glm( n_recaps / n_released ~ 1, 
                  family=binomial, 
                  data=tmp.df, 
                  weights=n_released ) 
      fit.AIC <- AIC(fit)
      
      cat(paste("df= ", 1, ", conv= ", fit$converged, " bound= ", fit$boundary, " AIC= ", round(fit.AIC, 4), "\n"))
      
      #added chunk from CAMPR eff_model.R that catches if model fit doesn't work
      if(!fit$converged | fit$boundary){
        cat("Constant (intercept-only) logistic model for efficiency failed. Using 'ROM+1' estimator. ")
        obs.mean<-tmp.df%>%
          summarise(ROM_eff=sum(na.omit(n_recaps+1))/sum(na.omit(n_released)+1))
        cat(paste(t,"'ROM+1' efficiency= ", trap_mean_eff$ROM_eff, "\n"))
        
        df$efficiency <- obs.mean
        
        fits[[trap]] <- fit
        eff.type[[trap]] <- 3
      } else {
        ####################################################################
        #2nd: fit increasingly complex models with increasing dof
        ####################################################################
        #skip the quadratic df =2
        #df = 3 = cubic model (no internal knots)
        #df = 4 = cubic spline w/ 1 internal knot at median
        #df = 5 = cubic spline w/ 2 internal knots at 0.33 and 0.66 of range
        #etc. (subtract 3 from df to get number of internal knots)
        cur.df <- 3
        repeat{
          
          #create spline basis for all days in season
          cur.bspl <- splines::bs( df$batch_date[ind.inside],#set current b-splines
                                   df=cur.df )
          
          #subset to only days with efficiency trials
          tmp.bs <- cur.bspl[!is.na(df$efficiency[ind.inside]),] 
          
          #fill the glm model with spline terms
          cur.fit <- glm( efficiency ~ tmp.bs, 
                          family=binomial,
                          data=tmp.df, 
                          weights=tmp.df$n_released ,#when n_released is larger, the efficiency a estimate is more precise, so the glm puts higher value in those because we don't assume constant variance
                          na.action = na.exclude, 
                          singular.ok = FALSE)
          
          #calculate AIC
          cur.AIC <- AIC(cur.fit)
          cat(paste("df= ", cur.df, ", conv= ", cur.fit$converged, " bound= ", cur.fit$boundary, " AIC= ", round(cur.AIC, 4), "\n"))
          
          #check results to see if need to break repeat loop
          if( !cur.fit$converged | cur.fit$boundary | cur.df > max.df.spline | cur.AIC > (fit.AIC + 4) ){
            break
          } else {
            fit <- cur.fit
            fit.AIC <- cur.AIC
            bspl <- cur.bspl
            cur.df <- cur.df + 1
          }
        }
        
        cat("\nFinal Efficiency model for trap: ", trap, "\n")
        print(summary(fit, disp=sum(residuals(fit, type="pearson")^2)/fit$df.residual))
        
        #make a design matrix for ease in calculating predictions.
        if( length(coef(fit)) <= 1 ){
          pred <- matrix( coef(fit), sum(ind.inside), 1 )
          X <- matrix( 1, sum(ind.inside), 1)
        } else {
          X <- cbind( 1, bspl )
          pred <- X %*% coef(fit)
        }
        
        #save X, and the dates at which we predict, for bootstrapping.

        all.dts[[trap]] <- df$batch_date[ind.inside]   
        
        #standard logistic prediction equation.  
        #"Pred" is all efficiencies for dates between min and max of trials.
        pred <- 1 / (1 + exp(-pred))  
        
        #replace all days efficiency (within trial days) with predicted
        df$efficiency_imputed[ind.inside] <- pred
        
        #use the mean of predictions for all dates outside trial period
        mean.p<-mean(pred,na.rm=T)
        df$efficiency_imputed[!ind.inside]<-mean.p
        
        #note: I believe the following is a bug and have commented it out
        #it is to reduce the dataset after expanding dates, but the
        #expansion was only done in the enhanced efficiency model chunk
        #df <- df[df$batchDate >= min.date.p & df$batchDate <= max.date.p,]
        fits[[trap]] <- fit
        
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
      df$eff_imputed <- factor( !ind, levels=c(T,F), labels=c(TRUE, FALSE))
      df$efficiency<-ifelse(df$eff_imputed==TRUE,
                            df$efficiency_imputed,df$efficiency)
      #if original trial efficiency is 0, impute, this is not done in
      #campR because everything is imputed
      #we have to do this because eff of 0 makes inf catch
      df <- df %>%
        mutate(
          eff_imputed = ifelse(efficiency == 0, TRUE, 
                               as.character(eff_imputed)),
          efficiency = ifelse(efficiency == 0, efficiency_imputed, efficiency)
        )
    }else{
      df$eff_imputed<-TRUE
      df$efficiency<-df$efficiency_imputed
    }
    df$trap_ID_decimal <- trap
    df$eff_imputed<-as.logical(df$eff_imputed)
    
    if(!is.na(fits[trap])){
      #say if we used enhanced efficiency model (in this chunk we don't)
      df$enhanced.eff <- rep("No",nrow(df))
      df <- df[,c("trap_ID","batch_date","n_released","n_recaps","eff_imputed","efficiency","enhanced.eff","trap_ID_decimal")]
    }
    
    if(is.na(fits[trap])){  
      #save the fit for bootstrapping.
      ans$models[[trap]] <- NA
      ans$all.X[[trap]] <- NA
      ans$all.ind.inside=all.ind.inside #provide first and last date per trap
      ans$all.dts[[trap]]<-df$batch_date[ind.inside] #provide dates that X predicts for
      ans$obs.data[[trap]]<-tmp.df
      ans$eff.type<-eff.type
      df<-df[,c("trap_ID","trap_ID_decimal","batch_date","n_released","n_recaps","eff_imputed")]
      df<-df%>%mutate(eff_imputed=NA,
                      efficiency=NA,
                      enhanced.eff=NA)
      ans$results<-ans$results%>%rbind(df)
    } else {
      #save the fit for bootstrapping.
      ans$models[[trap]] <- fit
      ans$all.X[[trap]] <- X
      ans$all.ind.inside=all.ind.inside #provide first and last date per trap
      ans$all.dts[[trap]]<-df$batch_date[ind.inside] #provide dates that X predicts for
      ans$obs.data[[trap]]<-tmp.df
      ans$eff.type<-eff.type
      
      ans$results<-ans$results%>%rbind(df)
      }
    
  }
  
  return(ans)
}
