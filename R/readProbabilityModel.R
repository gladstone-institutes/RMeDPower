readProbabilityModel <- function(jsonfile) {
  Prob_Model_data <- jsonlite::fromJSON(jsonfile)

  prob_model <- new("ProbabilityModel")

  prob_model@error_is_non_normal = Prob_Model_data$error_is_non_normal
  if(!is.null(Prob_Model_data$family_p)) prob_model@family_p = Prob_Model_data$family_p

  return(prob_model)
}
