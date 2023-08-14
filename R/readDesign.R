
readDesign <- function(jsonfile, data) {
  design_data <- jsonlite::fromJSON(jsonfile)

  design <- new("RMeDesign")

  design@response_column = design_data$response_column
  if(!is.null(design_data$covariate)) design@covariate = design_data$covariate
  design@condition_column = design_data$condition_column
  design@condition_is_categorical = design_data$condition_is_categorical
  design@experimental_columns = design_data$experimental_columns
  if(!is.null(design_data$crossed_columns)) design@crossed_columns = design_data$crossed_columns
  if(!is.null(design_data$total_column)) design@total_column = design_data$total_column
  design@outlier_alpha = design_data$outlier_alpha

  return(design)
}
