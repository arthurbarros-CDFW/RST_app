#' @title model_catch
#'
#' @description Fits Poisson GLM to catch data and imputes catch for missing periods. Called by est_catch().
#' Takes the following steps:
#' 1) Fits a null model to catch data for each trap sampling segment.
#' 2) Calculates AIC for model fitting.
#' 3) If enough real catch data available, increase GLM spline degrees of freedom until model fit is good enough
#'  (ie model doesn’t converge or AIC increases etc.).
#' 4) Once model is fit, impute catch by looping over batch date sampling gaps.
#' 
#' @param catch.df Input of catch data frame with periods of missing catch data.
#' 
#' @return ans List containing catch model fit results and catch data with imputed catch values.
#' Also contains X.miss, model matrix with intercept and predictor terms for generating catch predictions.

model_catch<-function(catch.df){
  
  traps<-unique(catch.df$trap_ID_decimal)
  
  ans<-list()
  ans$results<-NULL
  ans$models<-list()
  
  for(trap in traps){
    d<-catch.df%>%filter(trap_ID_decimal==trap)
    
    #before I tried filtering by na start_time, but this didn't
    #work because the "decimal site ID" traps after gaps still had
    #a start_time, and visitType doesn't work because sometimes
    #trapping starts after short gaps
    #solution: just drop the first record in each new ID_decimal
    #for the first batch date
    #also consider the "is_long_gap" field from the 
    #calculate_sampling_time()
    d<-d %>%
      slice(-which.min(batch_date))
    
    #sort data
    d <- d[ order(d$trap_ID_decimal, d$end_time), ]
    
    #filter out periods not sampling
    df<-d%>%filter(trap_status!="Not fishing")
    
    #define an offset in terms of hours.
    log.sample_hrs <-  log(as.numeric( df$total_sample_minutes/60 ))
    
    #fit null model, the gap catches are NA
    fit<-glm(total_catch~offset(log.sample_hrs),
             family=poisson,
             data=df)
    
    #Ripped from CAMPR:
    #If we allow decimal fish to enter the Poisson catch model, we cannot readily
    #calculate AIC, which requires a likelihood.  This is because the AIC function
    #in the Poisson family in R utilizes the gamma function to calculate log(n!),
    #i.e., the cumulant.  This function breaks down if n is not an integer.  To 
    #remedy this, explicitly calculate the likelihood in the case when 
    #unassd.sig.digit is a positive integer other than zero. 
    #Nemes, Gergő (2010), "New asymptotic expansion for the Gamma function", 
    #Archiv der Mathematik, 95 (2): 161–169.
    nemes2007 <- function(vec){
      ans <- 0.5*(log(2*pi) - log(vec+1)) + (vec+1)*(log(vec+1+(1/ ( 12*(vec+1) - (1/(10*(vec+1))) ) ))-1)
      return(ans)
    }
    
    #get the AIC for the null model, depending on how we want the decimals to work out. 
    #Arthur note: not sure why we are worried about decimals here, revisit
    if(unassd.sig.digit == 0){
      fit.AIC <- AIC(fit)
    } else {
      poissonCumulant <- nemes2007(fit$y)
      loglikelihood <- sum(fit$y*fit$linear.predictors - fit$fitted.values - poissonCumulant)
      fit.AIC <- -2*loglikelihood + 2*length(coefficients(fit))
      fit$aic <- fit.AIC
    }
    
    #assemble the beta coefficients to output. 
    coefs <- data.frame(t(rep(NA,50)))
    coefs[1,1:length(coefficients(fit))] <- coefficients(fit)
    colnames(coefs) <- seq(1,50,1) 
    
    #get number of days with catch data
    nGoodData <- length(!is.na(df$total_catch))  
    catch.fit.all <- cbind(data.frame(trap_ID_decimal=df$trap_ID_decimal[1], df=0, conv=fit$converged, bound=fit$boundary, AIC=round(fit.AIC,4), nGoodData=nGoodData),coefs)
    
    #fit the glm model, increasing df until model breaks
    cat(paste("Number of non-zero catches : ", sum(!is.na(df$total_catch) & (df$total_catch > 0)), "\n"))
    cat("Catch model fitting:\n")
    cat(paste("df= ", 0, ", conv= ", fit$converged, " bound= ", fit$boundary, " AIC= ", round(fit.AIC, 4), "\n"))
    
    #estimate a spline of increasing complexity like we did with the efficiency model
    #basically each iteration increase degrees of freedom in the spline
    #then check to see if model should break or continue
    
    #first check to ensure we have at least 10 good trapping visits
    if(sum(nGoodData & df$total_catch>0)>10){
      cur.df <- 1  
      
      repeat{
        #check to make sure enough data points for current complexity
        valid <- (cur.df <= floor(nGoodData / knotMesh))
        
        #ensure at least 15 non-NA trapping end times
        if(cur.df==1 & valid==TRUE){
          #pull julian data for linear spline term
          j<-as.numeric(format(df$end_time, "%j"))
          bs.sEnd<-matrix(j,ncol=1)
        }else if(cur.df==2 & valid==TRUE){
          j<-as.numeric(format(df$end_time, "%j"))
          #add quadratic spline term
          bs.sEnd<-cbind(Lin=j,Quad=j*j)
        }else if(cur.df>2 & valid==TRUE){
          #increase degrees of freedom in spline
          bs.sEnd<-splines::bs(df$end_time,df=cur.df)
        }
        
        #store new current fit in memory
        #model is essentially: total_catch ~ sampling time + julian date
        cur.fit <- tryCatch(glm( total_catch ~ offset(log.sample_hrs) + bs.sEnd,
                                 family=poisson, data=df ),
                            error=function(e) e )
        
        #overwrite aic with nemes2007 correction
        if(unassd.sig.digit>0){
          poissonCumulant<-nemes2007(cur.fit$y)
          loglikelihood<-sum(cur.fit$y*cur.fit$linear.predictors - cur.fit$fitted.values)-sum(poissonCumulant)
          cur.fit.AIC<- -2*loglikelihood + 2*length(coefficients(cur.fit))
          cur.fit$aic<-cur.fit.AIC
        }
        
        #ensure current fit is interpretable
        if(class(cur.fit)[1]=="simpleError"){
          cur.AIC=NA
          cat(paste0("Model fell apart;  tryCatch caught the output.  You may wish to investigate for trap ",trap,".\n"))
        }else{
          cur.AIC <- AIC( cur.fit )
          cat(paste("df= ", cur.df, ", conv= ", cur.fit$converged, " bound= ", cur.fit$boundary, " AIC= ", round(cur.AIC, 4), "\n"))
        }
        
        #compare new fit with previous fit, do we keep going or stop?
        if(is.na(cur.AIC)){
          break
        }else if( !cur.fit$converged | cur.fit$boundary | cur.df>15 | cur.AIC > (fit.AIC-2)){
          break
        }else{
          
          #assemble beta coefficients
          coefs<-data.frame(t(rep(NA,50)))
          coefs[1,1:length(coefficients((cur.fit)))]<-coefficients(cur.fit)
          colnames(coefs)<-seq(1,50,1)
          
          catch.fit<-cbind(data.frame(trap_ID_decimal=trap,
                                      df=cur.df,
                                      conv=cur.fit$converged,
                                      bound=cur.fit$boundary,
                                      AIC=round(cur.AIC,4),
                                      nGoodData=nGoodData),
                           coefs)
          
          fit<-cur.fit
          fit.AIC<-cur.AIC
          bs.end_time<-bs.sEnd
          cur.df<-cur.df+1 #increase df
        }
        #bind new fit to catch.fit.all
        catch.fit.all <- rbind(catch.fit.all,catch.fit)
      }
    }
    
    #done fitting model
    print(summary(fit, disp=sum(residuals(fit, type="pearson")^2)/fit$df.residual))
    
    #Next impute catch by looping over gaps one at a time
    #predict catch for max 24 hour periods
    
    d$catch_imputed <- FALSE
    
    #number of columns in smoother part excluding intercept
    degree<-length(coef(fit))-1
    if(degree<=2){
      nots<- -1
      b.knots<- -1
    } else {
      nots<- attr(bs.end_time, "knots")
      b.knots<-attr(bs.end_time, "Boundary.knots")
    }
    
    #setup for creation
    i<-1
    all.new.dat <- NULL
    all.gaplens <- NULL
    all.bdates <- NULL
    jason.new <- NULL
    
    #sweep through the df imputing catch data where needed
    repeat{
      if(i>=nrow(d)) break
      if(d$trap_status[i]=="fishing"){
        i<-i+1
      } else if(d$total_sample_minutes[i]<=(60*max.ok.gap)){
        #deal with small gaps (just remove)
        if(i>1){
          df$end_time[i-1]<-df$end_time[i]
          df$batch_date[i-1]<-df$batch_date[i]
          df<-df[-i,]
        }
      } else {
        
        #here we have missing values with period of time > max.ok.gap with no fish estimates
        
        #break up not fishing period into 24 hour chunks and remainder
        i.gapLens<-c(rep(24,
                         floor(d$total_sample_minutes[i]/(24*60))),
                     (d$total_sample_minutes[i]/60) %% 24)
        
        #make sure last remainder interval is >max.ok.gap
        ng<-length(i.gapLens)
        if(i.gapLens[ng]<=max.ok.gap){
          #add last small gap to prior interval
          i.gapLens[ng-1]<-i.gapLens[ng-1]+i.gapLens[ng]
          i.gapLens<-i.gapLens[-ng]
        }
        ng<-length(i.gapLens)
        
        #calculate the length of time for which we need an estimate
        #this creates a new chunk for each day in the gap
        sEnd<-d$start_time[i]+cumsum(i.gapLens*60*60)
        class(sEnd) <- class(d$start_time)
        attr(sEnd, "tzone") <- time_zone
        
        sStart<-d$start_time[i] +cumsum(c(0,i.gapLens[-ng])*60*60)
        class(sStart) <- class(sEnd)
        attr(sStart, "tzone") <- time_zone
        
        #the smoother
        if(degree<=2){
          bs.sEnd<- -1
        }else{
          bs.sEnd<-splines::bs(sEnd,knots=nots,
                               Boundary.knots = b.knots)
          dimnames(bs.sEnd)[[2]]<-paste("bs.end_time",
                                        dimnames(bs.sEnd)[[2]],
                                        sep="")
        }
        
        if(degree ==0){
          #mean model (intercept only)
          new.dat<-matrix(1, length(sEnd),1)
        }else if(degree==1){
          #linear model (intercept + julian day)
          j<-as.numeric(format(sEnd,"%j"))
          new.dat<-cbind(rep(1,length(sEnd)),j)
        }else if(degree==2){
          #quadratic model (intercept + julian day +julian day^2)
          j <- as.numeric( format(sEnd, "%j") )
          new.dat<-cbind( rep(1, length(sEnd)), j, j*j )
        } else {
          #cubic spline model (uses B-spline matrix computed earlier)
          new.dat <- cbind( 1, bs.sEnd )
        }
        if(degree > 2){
          bs.object <- list(
            knots = attr(bs.end_time, "knots"),
            Boundary.knots = attr(bs.end_time, "Boundary.knots"),
            degree = 3  # cubic splines
          )
          # Save this with the model
          fit$bs.object <- bs.object
        }
        
        #get the catch estimate based on model
        #parameters and length of time
        pred <- (new.dat %*% coef(fit)) + log(i.gapLens)
        pred <- exp(pred)
        
        #put things into blank database for output
        #here Trent uses a "trick" to build a data frame
        #with the same structure as the original
        new<-d[1:ng,]
        
        #put in imputed catch data
        new$total_catch<-as.numeric(pred)
        
        new$total_sample_minutes <- i.gapLens * 60
        new$end_time <- sEnd
        new$start_time <- sStart
        new$catch_imputed <- TRUE
        new$trap_status<-"Not fishing"
        new$trap_ID <- d$trap_ID[i]
        new$trap_ID_decimal <- d$trap_ID_decimal[i]
        new <- assign_batch_date(new , "end_time")
        
        #insert new data frame into d
        d<-rbind(d[1:(i-1),], new, d[(i+1):nrow(d),])
        
        #save imputed numbers for checking daily counts in baseTable.
        jason.new <- rbind(jason.new,new)
        
        #save model matrix used to make predictions for future bootstrapping
        all.new.dat<-rbind( all.new.dat, new.dat )
        all.gaplens <- c(all.gaplens, i.gapLens)
        
        all.bdates <- c(all.bdates, new$batch_date)  
        i <- i + ng + 1
      }
      #no break needed here
    }
    
    if( is.null( all.new.dat )){
      all.new.dat <- NA
    }
    
    ans$models[[trap]] <- fit
    
    ans$results<-ans$results%>%rbind(d)
    
    ans$X.miss[[trap]]<-all.new.dat
    
    ans$gaps[[trap]]=all.gaplens
    
    ans$batchDate.for.missings[[trap]]=(all.bdates)
  }
  return(ans)
}
