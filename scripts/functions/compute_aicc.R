#Burnham, K. P., & Anderson, D. R. (2002). Model selection and multimodel inference.
compute_aicc <- function(model) {
  n <- nobs(model)  # number of observations used in fitting
  k <- length(coef(model))  # number of parameters
  aic <- AIC(model)
  
  #AICc formula: AIC + (2*k*(k+1))/(n - k - 1)
  #only apply when n - k - 1 > 0
  if((n - k - 1) > 0) {
    aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
    return(aicc)
  } else {
    #if too many parameters for sample size, return Inf (model is invalid)
    return(Inf)
  }
}
