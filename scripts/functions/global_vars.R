#' @title global_vars
#'
#' @description Defines default variables used across all functions.
#' 
#' @param gap_threshold_days Sequential days of a sampling gap that define a new trapping segment (default: 7).
#' 
#' @param time_zone 	Timezone for date handling (default: America/Los_Angeles).
#' 
#' @param unassd.sig.digit Decimal handling for AIC correction in model selection (1 = use Nemes 2007 correction).
#' 
#' @param knotMesh Minimum data points per spline degree of freedom (default: 15).
#' 
#' @param max.ok.gap Maximum gap (hours) that can be ignored without imputation (default: 2).
#' 
#' @param bootstrap.CI.fx Function for bootstrap confidence intervals (default: "f.ci").
#' 
#' @param sum.by Time grouping variable to summarize passage data by (default: "week").


gap_threshold_days=7
time_zone<-"America/Los_Angeles"
unassd.sig.digit = 1
knotMesh = 15 #used to determine complexity of spline in model_catch.R
max.ok.gap=2 #maximum gap in sampling hours that is "okay"
bootstrap.CI.fx = "f.ci"
sum.by="week"