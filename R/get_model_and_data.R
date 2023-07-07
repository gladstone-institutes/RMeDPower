#' @title calculate_lmer_estimates
#'
#'
#' @description This function performs a linear mixed model analysis using lmer.
#'
#' Note: The current version does not accept categorical response variables, sample size parameters smaller than the observed samples size
#'
#' @import multtest
#' @import simr
#' @import lme4
#' @import lmerTest
#' @import readxl
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of variables related to experimental design such as "experiment", "plate", and "cell_line". "experiment" should come always first
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param total_column Set this column only when family_p="binomial" and it is equal to the total number of observations (number of cases plus number of controls) for a given number of cases
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param crossed_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param error_is_non_normal  error distribution of response . TRUE: Non-normal distribution , FALSE: Normal distribution (default).
#' @param family_p The type of distribution family to specify when the response is categorical. If family is "binary" then binary(link="log") is used, if family is "poisson" then poisson(link="logit") is used, if family is "poisson_log" then poisson(link=") log") is used.
#' @param na.action "complete": missing data is not allowed in all columns (default), "unique": missing data is not allowed only in condition, experimental, and response columns. Selecting "complete" removes an entire row when there is one or more missing values, which may affect the distribution of other features.
#'
#' @return A list of the linear mixed model result, original data, experimental column names, and residual values
#'
#' @export
#' @examples



get_model_and_data <- function(data, condition_column, experimental_columns, response_column, total_column = NULL, condition_is_categorical,
                                     crossed_columns=NA, error_is_non_normal=FALSE, family_p=NULL, na.action="complete"){



  ######input error handler
  if(!condition_column%in%colnames(data)){ print("condition_column should be one of the column names");return(NULL) }
  if(sum(experimental_columns%in%colnames(data))!=length(experimental_columns) ){ print("experimental_columns must match column names");return(NULL) }
  if(!response_column%in%colnames(data)){  print("response_column should be one of the column names");return(NULL) }
  if(is.null(condition_is_categorical) | !condition_is_categorical%in%c(TRUE,FALSE)){ print("condition_is_categorical must be TRUE or FALSE");return(NULL) }
  if(!is.na(crossed_columns)){if(sum(crossed_columns%in%colnames(data))!=length(crossed_columns) ){ print("crossed_columns must match column names");return(NULL) }}

  if(error_is_non_normal==TRUE){
    if(family_p != "negative_binomial")
      family_p=switch(family_p, "poisson" = poisson(link="log"), "binomial" = binomial(link="logit"), "bionomial_log" = binomial(link="log") )
    else
      family_p = list(family = "negative_binomial")
  }


  if(na.action=="complete"){

    notNAindex=which( rowSums(is.na(data)) == 0 )

  }else if(na.action=="unique"){

    notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column)])) == 0 )

  }



  fixed_global_variable_data=data[notNAindex,]


  cat("\n")
  print("__________________________________________________________________Summary of data:")
  print(summary(fixed_global_variable_data))
  cat("\n")

  colnames_original=colnames(fixed_global_variable_data)
  experimental_columns_index=NULL
  ####### assign categorical variables
  if(condition_is_categorical==TRUE) fixed_global_variable_data[,condition_column]=as.factor(fixed_global_variable_data[,condition_column])

  cat("\n")


  noncrossed_columns=NULL

  for(i in 1:length(experimental_columns)){
    fixed_global_variable_data[,experimental_columns[i]]=as.factor(fixed_global_variable_data[,experimental_columns[i]])
    experimental_columns_index=c(experimental_columns_index,which(colnames(fixed_global_variable_data)==experimental_columns[i]))
    colnames(fixed_global_variable_data)[experimental_columns_index[i]]=paste("experimental_column",i,sep="")

    if(i!=1&&!experimental_columns[i]%in%crossed_columns){
      noncrossed_columns=c(noncrossed_columns, paste("experimental_column",i,sep=""))
    }


    cat("\n")
    print(paste("_________________________________",experimental_columns[i]," is assigned to experimental_column",i,sep=""))
    cat("\n")
  }



  if(length(experimental_columns)>=2){
      for(r in 2:length(experimental_columns)){
      if(colnames(fixed_global_variable_data)[experimental_columns_index[r]]%in%noncrossed_columns){
        fixed_global_variable_data[,experimental_columns_index[r]]=paste(fixed_global_variable_data[,experimental_columns_index[r-1]], fixed_global_variable_data[,experimental_columns_index[r]],sep="_")
      }
    }

  }



  colnames(fixed_global_variable_data)[which(colnames(fixed_global_variable_data)==condition_column)]="condition_column"
  colnames(fixed_global_variable_data)[which(colnames(fixed_global_variable_data)==response_column)]="response_column"

  if(!is.null(total_column))
    colnames(fixed_global_variable_data)[which(colnames(fixed_global_variable_data)==total_column)]="total_column"


  ####### run the formula

  if(error_is_non_normal==FALSE){
    if(length(experimental_columns)==1){
      lmerFit <- lmerTest::lmer(response_column ~ condition_column + (1 | experimental_column1), data=fixed_global_variable_data)
    }else if(length(experimental_columns)==2){
      lmerFit <- lmerTest::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2), data=fixed_global_variable_data)
    }else if(length(experimental_columns)==3){
      lmerFit <- lmerTest::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=fixed_global_variable_data)
    }else if(length(experimental_columns)==4){
      lmerFit <- lmerTest::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=fixed_global_variable_data)
    }else if(length(experimental_columns)==5){
      lmerFit <- lmerTest::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=fixed_global_variable_data)
    }
  }else if(family_p$family == "binomial"){
    if(length(experimental_columns)==1){
      lmerFit <- lme4::glmer(cbind(response_column, (total_column - response_column)) ~ condition_column + (1 | experimental_column1), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==2){
      lmerFit <- lme4::glmer(cbind(response_column, (total_column - response_column)) ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==3){
      lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==4){
      lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==5){
      lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=fixed_global_variable_data, family=family_p)
    }
  }else if(family_p$family == "negative_binomial" & !is.null(total_column)){
    if(length(experimental_columns)==1){
      lmerFit <- lme4::glmer.nb(response_column ~ condition_column + (1 | experimental_column1) + offset(log(total_column)), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==2){
      lmerFit <- lme4::glmer.nb(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + offset(log(total_column)), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==3){
      lmerFit <- lme4::glmer.nb(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + offset(log(total_column)) , data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==4){
      lmerFit <- lme4::glmer.nb(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + offset(log(total_column)), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==5){
      lmerFit <- lme4::glmer.nb(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5) + offset(log(total_column)) , data=fixed_global_variable_data, family=family_p)
    }
  }else{
    if(length(experimental_columns)==1){
      lmerFit <- lme4::glmer(response_column ~ condition_column + (1 | experimental_column1), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==2){
      lmerFit <- lme4::glmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==3){
      lmerFit <- lme4::glmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==4){
      lmerFit <- lme4::glmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=fixed_global_variable_data, family=family_p)
    }else if(length(experimental_columns)==5){
      lmerFit <- lme4::glmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=fixed_global_variable_data, family=family_p)
    }
  }







  slmerFit <- summary(lmerFit)

  cat("\n")
  print("__________________________________________________________________Model statistics:")
  print(slmerFit)
  cat("\n")



  return(list(lmerFit, fixed_global_variable_data, colnames(fixed_global_variable_data)[experimental_columns_index], slmerFit$residuals))
}

