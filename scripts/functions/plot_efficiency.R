#' @title plot_efficiency
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

plot_efficiency<-function(impute_eff_results,
                          use_discharge,
                          impute_all){
  
  d<-impute_eff_results
  #filter out NA discharge
  if(use_discharge){
    discharge_data <- d %>%
      filter(!is.na(discharge) & !is.nan(discharge))
    scale_factor <- max(d$efficiency, na.rm = TRUE) / 
      max(d$discharge, na.rm = TRUE)
  }
  
  p_eff <- ggplot()
  
  #add conditional layers
  if (impute_all == FALSE) {
    p_eff <- p_eff +
      geom_point(data = d, 
                 aes(x = batch_date, y = efficiency_imputed, 
                     color = "Imputed Efficiency",
                     shape = "Imputed Efficiency")) +
      geom_point(data = d, 
                 aes(x = batch_date, y = unimputed_efficiency, 
                     color = "Unimputed Efficiency",
                     shape = "Unimputed Efficiency"))
  } else if (impute_all == TRUE) {
    p_eff <- p_eff +
      geom_point(data = d, 
                 aes(x = batch_date, y = efficiency_imputed, 
                     color = "Imputed Efficiency",
                     shape = "Imputed Efficiency"))
  }
  
  #add discharge points
  if (use_discharge == TRUE) {
    p_eff <- p_eff +
      geom_point(data = discharge_data, 
                 aes(x = batch_date, y = discharge * scale_factor, 
                     color = "Discharge",
                     shape = "Discharge"), 
                 alpha = 0.6, size = 1.5) +
      scale_color_manual(
        name="Values",
        values = c("Imputed Efficiency" = "#00BFC4", 
                   "Unimputed Efficiency" = "#F8766D",
                   "Discharge" = "black")
      ) +
      scale_shape_manual(
        name="Values",
        values = c("Imputed Efficiency" = 15,    
                   "Unimputed Efficiency" = 16,  
                   "Discharge" = 21)
      ) +
      scale_y_continuous(
        name = "Efficiency",
        sec.axis = sec_axis(
          ~ . / scale_factor,
          name = "Discharge (cfs)"
        )
      )
  } else{
    p_eff <- p_eff +
      scale_color_manual(
        name="Values",
        values = c("Imputed Efficiency" = "#00BFC4", 
                   "Unimputed Efficiency" = "#F8766D")
      ) +
      scale_shape_manual(
        name="Values",
        values = c("Imputed Efficiency" = 15,    
                   "Unimputed Efficiency" = 16)
      ) +
      scale_y_continuous(
        name = "Efficiency"
      )
  }
  
  p_eff <- p_eff +
    theme_bw() +
    facet_wrap(. ~ trap_ID_decimal)
  
  return(p_eff)
}

