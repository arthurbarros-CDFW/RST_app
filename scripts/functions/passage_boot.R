###############################
#function to bootstrap passage estimates for uncertainty
###############################

#bootstrap variability comes from two sources:
#imputed catch and imputed efficiency
passage_boot<-function(passage_data,
                       sum.by,
                       catch.fits,
                       catch.X.miss, 
                       catch.gapLens,
                       catch.bDates.miss,
                       eff.fits,eff.X,
                       eff.ind.inside,
                       eff.X.dates,
                       eff.X.obs.data,
                       eff.type,
                       survey_start,
                       survey_end,
                       R=100,
                       conf=0.95,
                       ci=T){
  bootstrap.CI.fx <- get("bootstrap.CI.fx",envir=.GlobalEnv)
  
  #run summarize_passage() to averages over traps and sums by sum.by group
  n.orig <- summarize_passage(passage_data, sum.by )
  
  
  n.len <- nrow(n.orig)
  n.passage_data <- nrow(passage_data)
  
  if( !ci ){ #check if ci is called for (default=T) if not return empty fields
    na <- rep(NA, n.len)
    ans <- data.frame(l = na, u= na)
    correct <- NA
  } else {
    n.traps <- length(catch.fits)
    
    #matrices to hold bootstrap iterations
    catch.pred <- matrix( passage_data$total_catch, nrow=nrow(passage_data), ncol=R )      
    eff.pred <- matrix( passage_data$efficiency, nrow=nrow(passage_data), ncol=R )
    
    cat(paste("n.traps=", n.traps, "\n"))
    cat(paste("summed by", sum.by, "\n"))
    cat(paste("n.len=", n.len, "\n"))
    
    #main iteration loops over each trap
    for(trap in 1:n.traps){
      trapID <- names(catch.fits)[trap]
      trap.ind <- passage_data$trap_ID_decimal == trapID
      
      cat(paste("trap=", trapID, "\n" ))
      
      ####################
      #generate random imputations of catch
      ####################
      ind <- which(trapID == names(catch.fits))
      if( length(ind) > 0 & all(!is.na(catch.X.miss[[trapID]]))){
        c.fit <- catch.fits[[trapID]]
        X <- catch.X.miss[[trapID]]
        gaps <- catch.gapLens[[trapID]]
        bd.miss <- catch.bDates.miss[[trapID]]
        
        cat("in bootstrap_passage.r (hit return)...")
        
        if( all(!is.na(bd.miss)) & all(!is.null(bd.miss)) ){
          
          #have some gaps in the data
          
          #estimate overdispersion parameter
          disp<-over_dispersion(c.fit,family="poisson",type="pearson")
          
          #scale variance-covariance matrix scaled by dispersion
          sig <- disp * vcov( c.fit )  
          
          cat(paste("...Poisson over-dispersion in catch model for trap ", trapID, " = ", disp, "\n"))
          cat("in bootstrap_passage.r")
          
          #model coefficients
          beta <- coef( c.fit )
          
          #catch if number of non-zero catches <10 or no catch at all
          if( (length(beta) == 1) & (beta[1] < -10)){
            rbeta <- matrix( -10, nrow=R, ncol=1 )
          }else{
            #generate random covariates
            #R random realizations of the beta vector. rbeta is R X (n coef)
            rbeta <- mvtnorm::rmvnorm(n=R, mean=beta, sigma=sig, method="chol")
            
            #predict catches using random coefficients
            pred<-X%*%t(rbeta)+log(gaps)
            pred <- exp(pred) 
            
            pred_summed <- apply( pred, 2, function(x,bd){tapply(x,bd,sum)}, bd=bd.miss )
            pred_summed <- matrix( unlist(pred_summed), nrow=length(unique(bd.miss)), ncol=R )
            #the unique call above is making pred less than the number of imputes necessary, I think I fixed this.
            
            unique_bd <- unique(bd.miss)
            
            rows_per_date <- split(which(trap.ind & (passage_data$p.c.imputed > 0)), 
                                   bd.miss)
            
            pred_expanded <- matrix(NA, nrow=sum(sapply(rows_per_date, length)), ncol=R)
            
            # Fill in the expanded matrix
            row_counter <- 1
            for(i in 1:length(unique_bd)) {
              date <- unique_bd[i]
              rows <- rows_per_date[[as.character(date)]]
              n_rows <- length(rows)
              
              if(n_rows > 0) {
                # Divide the daily total equally among all rows for that date
                daily_total <- pred_summed[i, , drop = FALSE]  # Keep as matrix
                daily_per_row <- daily_total / n_rows
                
                # Fill in the expanded matrix
                end_row <- row_counter + n_rows - 1
                pred_expanded[row_counter:end_row, ] <- matrix(rep(daily_per_row, each = n_rows), 
                                                               nrow = n_rows, ncol = R)
                row_counter <- end_row + 1
              }
            }
            
            #now use pred_expanded instead of pred
            pred <- pred_expanded
            
            #make catch matrix the correct size by including the observed counts
            ind.mat <- matrix( trap.ind & (passage_data$p.c.imputed > 0), 
                               nrow=n.passage_data, ncol=R )
            
            #for ind.mat=T index rows, replace with predicted values
            catch.pred[ind.mat] <- pred
            
          }
        }
      }
        
        ####################
        #generate random imputations of efficiency
        ####################
        ind <- which( trapID == names(eff.fits) )
        if( length(ind) > 0 ){
          e.fit <- eff.fits[[trapID]]
          e.X <- eff.X[[trapID]]
          eff.obs.data <- eff.X.obs.data[[trapID]]
          e.type <- eff.type[[trapID]]
          
          #The e.X design matrix should have columns in order from model_eff.R
          #but here we can double check and reset order on e.X, which is easier 
          #to manipulate than order of e.fit
          #AB NOTE: this only seems to matter for the enhanced efficiency model with covariates
          
          if( (length(e.X) > 1) & !is.na(e.type) ){  #exclude cases when e.X == NA.  
            if(ncol(e.X) > 1 & e.type == 5){
              cat(paste0("Sorting variables ",colnames(e.X)," for e.X in bootstrap.\n"))
              timeVar <- sort(colnames(e.X)[grepl("time",colnames(e.X),fixed=TRUE)])  
              fit.Vars <- names(coef(e.fit))
              notTimeVar <- fit.Vars[!(grepl("tmp.bs",fit.Vars,fixed=TRUE))]
              notTimeVar <- notTimeVar[notTimeVar != "(Intercept)"]
              thisOrder <- c("Intercept",timeVar,notTimeVar)
              e.X <- e.X[,thisOrder]
            }
          }
          
          #another 2-vector of the first and last efficiency trials.
          e.ind <- eff.ind.inside[[trapID]]  
          
          #indicator for days inside the efficiency season.
          e.ind <- (e.ind[1] <= passage_data$batch_date) & (passage_data$batch_date <= e.ind[2]) & trap.ind 
          
          #vector of dates inside efficiency season
          e.dts <- eff.X.dates[[trapID]] 
          
          if( (!is.list(e.fit) & e.type != 5) | length(e.fit) == 0 | (e.type == 5 & length(e.fit) == 0) ){
            #in CAMPR this is an in progress error catch?
            print("no efficiency data found for bootstrapping")
          } else {
            
            #variance matrix: check if less than min_sample_size 
            if( e.fit$df.residual == 0 | nrow(e.fit$data) < min_sample_size ){
              disp <- 1
            } else {
              #estimate overdispersion parameter
              disp <- over_dispersion(model=e.fit,family="binomial",type="pearson")
            }
            
            cat(paste0("...Binomial over-dispersion in efficiency model for trap ",trapID," = ",disp, ".\n"))
            
            #scale variance-covariance matrix by overdispersion
            sig <- disp * vcov( e.fit )   
            
            #pull coefficients.
            beta <- coef( e.fit )
            
            #generate R random coefficients of the beta vector
            if( length(coef(e.fit)) == 1 ){ #interecept only model
              
              #if using enh efficiency model, not yet implemented
              if(e.type == 5){
                
                #bias adjustment based on the data that went into the enh eff estimation.  
                X <- rep(1,length(e.X))
                p <- (sum(e.fit$data$n_recaps) + 1) / (sum(e.fit$data$n_released) + 1)
                w <- e.fit$data$n_released*rep(p*(1 - p),length(X))
                
                diagonal <- matrix(rep(0,length(X)*length(X)),length(X),length(X))
                for(i in 1:length(X)){
                  diagonal[i,i] <- w[i]
                }
                
                #"matrix" is 1x1 here by design...only doing this for intercept-only models.
                sig <- as.matrix(disp*( solve(t(X) %*% diagonal %*% X) ))
                rbeta <- mvtnorm::rmvnorm(n=R, mean=log(p/(1-p)),sigma=sig,method="chol")
                
              }
              
              else {
                X <- rep(1,length(eff.obs.data$n_recaps))
                p <- (sum(eff.obs.data$n_recaps) + 1) / (sum(eff.obs.data$n_released) + 1) #ROM estimator for intercept only model
                w <- eff.obs.data$n_released*rep(p*(1 - p),length(X))
                
                diagonal <- matrix(rep(0,length(X)*length(X)),length(X),length(X))
                for(i in 1:length(X)){
                  diagonal[i,i] <- w[i]
                }
                
                #"matrix" is 1x1 here by design...only doing this for intercept-only models.
                sig <- as.matrix(disp*( solve(t(X) %*% diagonal %*% X) ))
                rbeta <- mvtnorm::rmvnorm(n=R, mean=log(p/(1-p)),sigma=sig,method="chol")
              }
              
            } else {#all other models non intercept only models way easier
              rbeta <- mvtnorm::rmvnorm(n=R, mean=beta, sigma=sig, method="chol") 
            }
            
            #predict efficiency using random coefficients
            pred<-(e.X %*% t(rbeta))
            pred <- 1 / (1 + exp(-pred))
            
            
            #use mean-predicted efficiency for times outside first and last trials.  
            ind.mat <- matrix( trap.ind, nrow=n.passage_data, ncol=R )
            e.means <- matrix( colMeans( pred ), byrow=T, nrow=sum(trap.ind), ncol=R )
            
            #this is complicated, but we have to line up the catch dates with 
            #the efficiency dates.  Because length of seasons vary, this is necessary.
            df.c <- data.frame(batch_date=format(passage_data$batch_date[trap.ind]), 
                               in.catch = TRUE, stringsAsFactors=FALSE )
            df.e <- data.frame(batch_date=format(e.dts), in.eff = TRUE, stringsAsFactors=FALSE )
            
            df.ce <- merge( df.c, df.e, all.x=TRUE )
            df.ec <- merge( df.e, df.c, all.x=TRUE )
            
            df.ec$in.catch[ is.na(df.ec$in.catch) ] <- FALSE
            df.ce$in.eff[ is.na(df.ce$in.eff) ] <- FALSE
            
            df.ce<-unique(df.ce)
            df.ec<-unique(df.ec)
            
            #predictions that are in the catch data set,
            pred <- pred[ df.ec$in.catch, ]   
            
            #indices where in.eff is TRUE
            eff_indices <- which(df.ce$in.eff) 
            
            if(length(eff_indices) != nrow(pred)) {
              cat(paste("Warning: eff_indices length =", length(eff_indices), 
                        "pred rows =", nrow(pred), "\n"))
              
              # try to align them - pred should correspond to the efficiency prediction dates
              # that are also in catch dates
              if(length(eff_dates) == nrow(pred)) {
                # Create a mapping
                eff_to_catch <- match(format(eff_dates), format(df.c$batch_date[df.ce$in.eff]))
                valid_matches <- !is.na(eff_to_catch)
                
                if(sum(valid_matches) > 0) {
                  # Assign only the matched rows
                  e.means[eff_indices[eff_to_catch[valid_matches]], ] <- pred[valid_matches, ]
                }
              }
            }else {
              #direct assignment works
              e.means[eff_indices, ] <- pred
            }
            
            #assign mean outside of eff.ind.inside[[trap]], and efficiency model inside season.
            eff.pred[ind.mat] <- e.means    
          }
      }
        
      cat("...BS complete\n")
    }
    
    pos <- 1
    envir <- as.environment(pos)
    assign("catch.pred", catch.pred, pos=envir)
    assign("eff.pred", eff.pred, pos=envir)
    
    #estimate passage in boot
    #matrices eff.pred and catch.pred are same size so just divide
    test <- ifelse(passage_data$p.c.imputed > 0 & passage_data$p.c.imputed < 1,
                   passage_data$total_catch - (passage_data$total_catch*passage_data$p.c.imputed),0)
    catch.pred <- apply(catch.pred,2,function(x) x + test)
    pass.pred <- catch.pred / eff.pred    #create replicates of passage prediction
    
    #   ---- Now, average over traps
    #   ---- At this point, pass.pred is a (n.batch.day*n.trap) X R matrix, with each cell containing the passage estimated
    #   ---- at a particular trap at the site for a particular batch day for a particular iteration of the bootstrap.
    #   ---- Row dimension of list items corresponds to (batch days x trap), columns correspond to iterations.
    #   ---- We now need to average the cells over the traps, and summarize by time.  Do this by calling 
    #   ---- summarize_passage on each column.
    
    #   ---- Internal function to summarize catch by s.by by applying
    #   ---- F.summarize to every column of pass.
    internal.sumize.pass <- function(p, s.by, bd){
      
      #   p <- c.pred
      df <- data.frame( batch_date=bd, passage=p, p.c.imputed=1 )
      
      #summarize passage for each rep
      n <- summarize_passage( df, s.by )
      n$passage
    }
    
    #apply above internal function to summarize passage for each boot iteration
    pass <- apply( pass.pred, 2, internal.sumize.pass,
                   s.by=sum.by, bd=passage_data$batch_date)
    pass <- matrix( unlist(pass), n.len, R )
    
    #another internal function to compute bias corrected boostrap CIs
    f.bias.acc.ci<-function(x,alpha,x.orig){
      p <- mean( x > rep(x.orig$passage,R), na.rm=TRUE)
      z.0 <- qnorm( 1 - p )
      z.alpha <- qnorm( 1 - (alpha/2))
      p.L <- pnorm( 2*z.0 - z.alpha )
      p.H <- pnorm( 2*z.0 + z.alpha )
      ci <- quantile( x[ !is.na(x) & (x < Inf) ], p=c(p.L, p.H) )
      ci
    }
    
    #"regular" bootstrap CIs based on defined CIs
    f.ci <- function( x, alpha, x.orig ){
      ci <- quantile( x[ !is.na(x) & (x < Inf) ], p=c(alpha/2, 1-(alpha/2)) )
      ci
    }
    
    #get (1 - alpha)% confidence bounds via GlobalVars bootstrap.CI.fx.
    #The below just checks to see what the default value for bootstrap.CI.fx
    #is in the globalvars. Default is "f.ci", so applies that default function
    #to calculate bootstrap CIs, we can change it to used bias corrected in future
    ans <- tryCatch({
      apply( pass, 1, get(bootstrap.CI.fx), alpha=(1-conf), x.orig=n.orig )
    },
    error=function(cond) {
      stop("I don't see a value for GlobalVars function bootstrap.CI.fx. Examine GlobalVars() for valid value(s).\n")
      # message(cond)
      return(NA)
    })
    
    ans <- as.data.frame(t(matrix( unlist(ans), 2, n.len )))
    
    correct <- apply(pass,1,function(x) sd(x))
    
  }
  names(ans) <- paste0( c("lower.", "upper."), conf*100 )
  
  ans <- data.frame( n.orig, ans, error=correct, stringsAsFactors=F )
  return(ans)
}
