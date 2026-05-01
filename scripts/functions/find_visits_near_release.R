###############################
#function finds sampling visits within 36 hours
#of a efficiency trial release
###############################
find_visits_near_release <- function(release_time, visit_df) {
  visit_df %>%
    filter(visit_datetime >= release_time & 
             visit_datetime <= release_time + hours(36))
}