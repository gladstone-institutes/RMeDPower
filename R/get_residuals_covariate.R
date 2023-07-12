#' @title get_residuals_covariate
#'
#' @description This function retrieve residual values from an lmerFit summary object and plot residual values by condition_column
#'
#' Note: The current version does not accept categorical response variables, sample size parameters smaller than the observed samples size
#'
#' @import multtest
#' @import lme4
#' @import lmerTest
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of variables related to experimental design such as "experiment", "plate", and "cell_line". They should be in order, for example, "experiment" should always come first .
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param total_column Set this column only when family_p="binomial" and it is equal to the total number of observations (number of cases plus number of controls) for a given number of cases
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param covariate The name of the covariate to control in the regression model
#' @param crossed_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param error_is_non_normal Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param family_p The type of distribution family to specify when the response is categorical. If family is "binary" then binary(link="log") is used, if family is "poisson" then poisson(link="logit") is used, if family is "poisson_log" then poisson(link=") log") is used.
#' @param na.action "complete": missing data is not allowed in all columns (default), "unique": missing data is not allowed only in condition, experimental, and response columns. Selecting "complete" removes an entire row when there is one or more missing values, which may affect the distribution of other features.
#' @return The adjusted residual values for the experimental variables will be returned along with the original input data. A boxplot of residual values and a plot of 95% confidence interval mean residual values adjusted for the experimental variables will be returned.
#'
#' @export
#'
#' @examples result=get_residuals_covariate(data=RMeDPower_data1,
#' @examples condition_column="classification",
#' @examples experimental_columns=c("experiment", "line"),
#' @examples response_column="cell_size1",
#' @examples condition_is_categorical=TRUE,
#' @examples covariate="covariate",
#' @examples crossed_columns = "line",
#' @examples error_is_non_normal=FALSE)


get_residuals_covariate<-function(data, condition_column, experimental_columns, response_column, total_column, condition_is_categorical, covariate=NA,
                        crossed_columns=NA, error_is_non_normal=FALSE, family_p=NULL, na.action="complete"){



  ######input error handler
  if(!is.na(covariate) & !covariate%in%colnames(data)){ print("covariate should be NA or one of the column names");return(NULL) }

  if(!condition_column%in%colnames(data)){ print("condition_column should be one of the column names");return(NULL) }
  if(sum(experimental_columns%in%colnames(data))!=length(experimental_columns) ){ print("experimental_columns must match column names");return(NULL) }

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

    if(is.na(covariate)) notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column, covariate)])) == 0 )
    else notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column)])) == 0 )

  }

  Data=data[notNAindex,]




  cat("\n")
  print("__________________________________________________________________Summary of data:")
  print(summary(Data))
  cat("\n")

  colnames_original=colnames(Data)
  experimental_columns_index=NULL
  ####### assign categorical variables
  if(condition_is_categorical==TRUE) Data[,condition_column]=as.factor(Data[,condition_column])

  cat("\n")


  noncrossed_columns=NULL

  for(i in 1:length(experimental_columns)){
    Data[,experimental_columns[i]]=as.factor(Data[,experimental_columns[i]])
    experimental_columns_index=c(experimental_columns_index,which(colnames(Data)==experimental_columns[i]))
    colnames(Data)[experimental_columns_index[i]]=paste("experimental_column",i,sep="")

    if(i!=1&&!experimental_columns[i]%in%crossed_columns){
      noncrossed_columns=c(noncrossed_columns, paste("experimental_column",i,sep=""))
    }


    cat("\n")
    print(paste("_________________________________",experimental_columns[i]," is assigned to experimental_column",i,sep=""))
    cat("\n")
  }


  if(length(experimental_columns)>=2){
        for(r in 2:length(experimental_columns)){
      if(colnames(Data)[experimental_columns_index[r]]%in%noncrossed_columns){
        Data[,experimental_columns_index[r]]=paste(Data[,experimental_columns_index[r-1]], Data[,experimental_columns_index[r]],sep="_")
      }
    }
  }




  colnames(Data)[which(colnames(Data)==condition_column)]="condition_column"
  colnames(Data)[which(colnames(Data)==response_column)]="response_column"
  if(!is.na(covariate)) colnames(Data)[which(colnames(Data)==covariate)]="covariate"

  if(!is.null(total_column))
    colnames(Data)[which(colnames(Data)==total_column)]="total_column"

  ####### run the formula

  if(is.na(covariate)){
    if(error_is_non_normal==FALSE){
      if(length(experimental_columns)==1){
        lmerFit <- lmerTest::lmer(response_column ~  (1 | experimental_column1), data=Data)
      }else if(length(experimental_columns)==2){
        lmerFit <- lmerTest::lmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2), data=Data)
      }else if(length(experimental_columns)==3){
        lmerFit <- lmerTest::lmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data)
      }else if(length(experimental_columns)==4){
        lmerFit <- lmerTest::lmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data)
      }else if(length(experimental_columns)==5){
        lmerFit <- lmerTest::lmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data)
      }
    }else if(family_p$family == "binomial"){
        if(length(experimental_columns)==1){
          lmerFit <- lme4::glmer(cbind(response_column, (total_column - response_column)) ~ (1 | experimental_column1), data=Data, family=family_p)
        }else if(length(experimental_columns)==2){
          lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ (1 | experimental_column1) + (1 | experimental_column2), data=Data, family=family_p)
        }else if(length(experimental_columns)==3){
          lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data, family=family_p)
        }else if(length(experimental_columns)==4){
          lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data, family=family_p)
        }else if(length(experimental_columns)==5){
          lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data, family=family_p)
        }
      }else if(family_p$family == "negative_binomial" & !is.null(total_column)){
        if(length(experimental_columns)==1){
          lmerFit <- lme4::glmer.nb(response_column ~  (1 | experimental_column1) + offset(log(total_column)), data=Data, family=family_p)
        }else if(length(experimental_columns)==2){
          lmerFit <- lme4::glmer.nb(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + offset(log(total_column)), data=Data, family=family_p)
        }else if(length(experimental_columns)==3){
          lmerFit <- lme4::glmer.nb(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + offset(log(total_column)) , data=Data, family=family_p)
        }else if(length(experimental_columns)==4){
          lmerFit <- lme4::glmer.nb(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + offset(log(total_column)), data=Data, family=family_p)
        }else if(length(experimental_columns)==5){
          lmerFit <- lme4::glmer.nb(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5) + offset(log(total_column)) , data=Data, family=family_p)
        }
      }else{
      if(length(experimental_columns)==1){
        lmerFit <- lme4::glmer(response_column ~  (1 | experimental_column1), data=Data, family=family_p)
      }else if(length(experimental_columns)==2){
        lmerFit <- lme4::glmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2), data=Data, family=family_p)
      }else if(length(experimental_columns)==3){
        lmerFit <- lme4::glmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data, family=family_p)
      }else if(length(experimental_columns)==4){
        lmerFit <- lme4::glmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data, family=family_p)
      }else if(length(experimental_columns)==5){
        lmerFit <- lme4::glmer(response_column ~  (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data, family=family_p)
      }
    }
  }else{
    if(error_is_non_normal==FALSE){
      if(length(experimental_columns)==1){
        lmerFit <- lmerTest::lmer(response_column ~  covariate + (1 | experimental_column1), data=Data)
      }else if(length(experimental_columns)==2){
        lmerFit <- lmerTest::lmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2), data=Data)
      }else if(length(experimental_columns)==3){
        lmerFit <- lmerTest::lmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data)
      }else if(length(experimental_columns)==4){
        lmerFit <- lmerTest::lmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data)
      }else if(length(experimental_columns)==5){
        lmerFit <- lmerTest::lmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data)
      }
    }else if(family_p$family == "binomial"){
      if(length(experimental_columns)==1){
        lmerFit <- lme4::glmer(cbind(response_column, (total_column - response_column)) ~ covariate + (1 | experimental_column1), data=Data, family=family_p)
      }else if(length(experimental_columns)==2){
        lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ covariate + (1 | experimental_column1) + (1 | experimental_column2), data=Data, family=family_p)
      }else if(length(experimental_columns)==3){
        lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data, family=family_p)
      }else if(length(experimental_columns)==4){
        lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data, family=family_p)
      }else if(length(experimental_columns)==5){
        lmerFit <- lme4::glmer(cbind(response_column, total_column - response_column) ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data, family=family_p)
      }
    }else if(family_p$family == "negative_binomial" & !is.null(total_column)){
      if(length(experimental_columns)==1){
        lmerFit <- lme4::glmer.nb(response_column ~ covariate + (1 | experimental_column1) + offset(log(total_column)), data=Data, family=family_p)
      }else if(length(experimental_columns)==2){
        lmerFit <- lme4::glmer.nb(response_column ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + offset(log(total_column)), data=Data, family=family_p)
      }else if(length(experimental_columns)==3){
        lmerFit <- lme4::glmer.nb(response_column ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + offset(log(total_column)) , data=Data, family=family_p)
      }else if(length(experimental_columns)==4){
        lmerFit <- lme4::glmer.nb(response_column ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + offset(log(total_column)), data=Data, family=family_p)
      }else if(length(experimental_columns)==5){
        lmerFit <- lme4::glmer.nb(response_column ~ covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5) + offset(log(total_column)) , data=Data, family=family_p)
      }
    }else{
      if(length(experimental_columns)==1){
        lmerFit <- lme4::glmer(response_column ~  covariate + (1 | experimental_column1), data=Data, family=family_p)
      }else if(length(experimental_columns)==2){
        lmerFit <- lme4::glmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2), data=Data, family=family_p)
      }else if(length(experimental_columns)==3){
        lmerFit <- lme4::glmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data, family=family_p)
      }else if(length(experimental_columns)==4){
        lmerFit <- lme4::glmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data, family=family_p)
      }else if(length(experimental_columns)==5){
        lmerFit <- lme4::glmer(response_column ~  covariate + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data, family=family_p)
      }
    }
  }







  slmerFit <- summary(lmerFit)
  cat("\n")
  print("__________________________________________________________________Model statistics:")
  print(slmerFit)
  cat("\n")

  residuals=slmerFit$residuals


  Data = cbind(residuals, Data)


  colnames(Data)[1] = "residual"
  #Data[,"condition_column"]=as.numeric(as.character(Data[,"condition_column"]))

  if(is.na(covariate)) {
    Data_sum <- Data %>%
      dplyr::group_by(condition_column, experimental_column1) %>%
      dplyr::summarise(mean_residual1 = mean(residual), med_residual1 = median(residual))
  }
  if(!is.na(covariate)){
    Data_sum <- Data %>%
      dplyr::group_by(covariate, condition_column, experimental_column1) %>%
      dplyr::summarise(mean_residual1 = mean(residual), med_residual1 = median(residual))
  }


  if(condition_is_categorical==TRUE){

    if(!is.na(covariate)) {
      gp=ggplot2::ggplot(Data_sum, aes(x=condition_column, y=med_residual1, color = covariate)) +
        geom_boxplot(position = position_dodge(width = 0.7), outlier.shape = NA) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.7)) +
        ylab("median residual") +
        theme_bw() +
        theme(axis.title   = element_text(face  = "bold"))
    }else{
      gp=ggplot2::ggplot(Data_sum, aes(x=condition_column, y=med_residual1)) +
        geom_boxplot(position = position_dodge(width = 0.5), outlier.shape = NA) +
        geom_point() +
        ylab("median residual") +
        theme_bw() +
        theme(axis.title   = element_text(face  = "bold"))
    }

    print(gp)
  }else{


    # plot(Data[,"condition_column"], Data[,"residual"], xlab=condition_column, ylab=paste0(response_column, " Residual Value"), main=NULL)
    # abline(lm( as.formula( paste0( "residual ~  condition_column") ), data=Data),col='blue')

    p <- ggplot(Data, aes(x=condition_column, y=residual, color=experimental_column1)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw()
    print(p)




    }

  return(Data)
}
