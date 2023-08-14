readPowerParams <- function(jsonfile) {
  Power_Param_data <- jsonlite::fromJSON(jsonfile)

  power_param <- new("PowerParams")

  power_param@target_columns = Power_Param_data$target_columns
  power_param@power_curve = Power_Param_data$power_curve
  power_param@nsimn = Power_Param_data$nsimn
  power_param@levels = Power_Param_data$levels
  power_param@alpha = Power_Param_data$alpha
  if(!is.null(power_param@max_size)) power_param@max_size = Power_Param_data$max_size
  if(!is.null(power_param@breaks)) power_param@breaks = Power_Param_data$breaks
  if(!is.null(power_param@effect_size)) power_param@effect_size = Power_Param_data$effect_size
  if(!is.null(power_param@icc)) power_param@icc = Power_Param_data$icc
  return(power_param)
}
