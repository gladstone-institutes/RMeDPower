readStatModel <- function(jsonfile) {
  Stat_Model_data <- jsonlite::fromJSON(jsonfile)

  stat_model <- new("StatModel")

  stat_model@error_is_non_normal = Stat_Model_data$error_is_non_normal
  if(!is.null(Stat_Model_data$family_p)) stat_model@family_p = Stat_Model_data$family_p

  return(stat_model)
}
