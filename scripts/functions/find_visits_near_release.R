#' @title find_visits_near_release
#'
#' @description Finds trap visits within 36 hours of a release. Called by est_efficiency() to map to individual release records.
#' 
#' @param release_time  Record of releases to be used for efficiency estimates.
#' 
#' @param visit_df For each release assigns visits from visit_df.
#' 
#' @return visit_df Returns back visit_df filtered to visits within 36 hours of release times.

find_visits_near_release <- function(release_time, visit_df) {
  visit_df %>%
    filter(visit_datetime >= release_time & 
             visit_datetime <= release_time + hours(36))
}