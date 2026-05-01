###############################
#function to get residuals of overdispersion
###############################

over_dispersion<-function(model,family,type){
  resids <- residuals(model,type)
  
  if(family == "binomial"){
    toss.ind <- (abs(resids) > 8)
    resids <- resids[!toss.ind]
    disp <- sum( resids*resids ) / (model$df.residual - sum(toss.ind))
    if( disp < 1.0 ){
      disp <- 1.0
    }
  }else if(family == "poisson"){
    qrds <- quantile( resids, p=c(.2, .8))
    toss.ind <- (resids < qrds[1]) | (qrds[2] < resids)
    resids <- resids[!toss.ind]
    disp <- sum( resids*resids ) / (model$df.residual - sum(toss.ind))
    if( disp < 1.0 ){
      disp <- 1.0
    }
  }
  return(disp)
  
}
