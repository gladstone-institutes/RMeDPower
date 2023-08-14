#' @title transform_data2_covariate
#'
#' @description This functions makes quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values after removing outliers. To detect outliers, the function uses Rosner's test.
#'
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of the variable related to experimental design such as "experiment", "plate", and "cell_line". They should be in order, for example, "experiment" should always come first .
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param total_column Set this column only when family_p="binomial" and it is equal to the total number of observations (number of cases plus number of controls) for a given number of cases
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param covariate The name of the covariate to control in the regression model
#' @param method The method used to detect outliers. "rosner" (default) runs Rosner's test and "cook" runs Cook's distance.
#' @param crossed_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param error_is_non_normal Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param family_p The type of distribution family to specify when the response is categorical. If family is "binary" then binary(link="log") is used, if family is "poisson" then poisson(link="logit") is used, if family is "poisson_log" then poisson(link=") log") is used.
#' @param alpha numeric scalar between 0 and 1 indicating the Type I error associated with the test of outliers
#' @param na.action "complete": missing data is not allowed in all columns (default), "unique": missing data is not allowed only in condition, experimental, and response columns. Selecting "complete" removes an entire row when there is one or more missing values, which may affect the distribution of other features.
#'
#' @return For continuous data, the function returns quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values. For discrete data, it returns a histogram. If "rosner" is chosen, a matrix with updated feature values after transformation will be returned. If "cook" is choose, a list with  a matrix with updated feature values after transformation will be returned,
#'
#' @export
#'
#' @examples result=transform_data2(data=data, condition_column="classif", experimental_columns=c("experiment","line"), response_column="feature", condition_is_categorical=TRUE, error_is_non_normal=FALSE, alpha=0.05, crossed_columns = "line", method="cook", na.action="complete")
#' @examples result=transform_data2(data=data, condition_column="classif", experimental_columns=c("experiment","line"), response_column="feature", condition_is_categorical=TRUE, error_is_non_normal=FALSE, alpha=0.05, crossed_columns = "line", method="cook", na.action="complete")
#' @examples result=transform_data2(data=data, condition_column="classif", experimental_columns=c("experiment","line"), response_column="feature", condition_is_categorical=TRUE, error_is_non_normal=TRUE, family_p="poisson", alpha=0.05, crossed_columns = "line", method="cook", na.action="complete")

diagnoseDataModel <- function(data, design, model) {
  transform_data2_covariate(data,
                  condition_column = design@condition_column,
                  experimental_columns = design@experimental_columns,
                  response_column = design@response_column,
                  total_column = design@total_column,
                  condition_is_categorical = design@condition_is_categorical,
                  covariate= design@covariate,
                  crossed_columns = design@crossed_columns,
                  error_is_non_normal=model@error_is_non_normal,
                  family_p=model@family_p,
                  alpha=design@outlier_alpha,
                  na.action=design@na_action)

}
