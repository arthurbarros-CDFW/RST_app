#' @title over_dispersion
#'
#' @description Calculate overdispersion for Poisson (catch) and Binomial (efficiency) models. Called by passage_boot().
#' 
#' @param model model fits developed from model_catch() and model_efficiency().
#' 
#' @param family Poisson for catch models, Binomial for efficiency models.
#' 
#' @param type Pearson residual type, measuring distance between predicted and actual values of model.
#' 
#' @return disp Returns dispersion metrics from models.

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
