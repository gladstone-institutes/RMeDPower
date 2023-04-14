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
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param covariate The name of the covariate to control in the regression model
#' @param repeatable_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param response_is_categorical Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
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
#' @examples repeatable_columns = "line",
#' @examples response_is_categorical=FALSE)


get_residuals_covariate<-function(data, condition_column, experimental_columns, response_column, condition_is_categorical, covariate=NA,
                        repeatable_columns=NA, response_is_categorical=FALSE, family_p=NULL, na.action="complete"){



  ######input error handler
  if(!is.na(covariate) & !covariate%in%colnames(data)){ print("covariate should be NA or one of the column names");return(NULL) }

  if(!condition_column%in%colnames(data)){ print("condition_column should be one of the column names");return(NULL) }
  if(sum(experimental_columns%in%colnames(data))!=length(experimental_columns) ){ print("experimental_columns must match column names");return(NULL) }

  if(is.null(condition_is_categorical) | !condition_is_categorical%in%c(TRUE,FALSE)){ print("condition_is_categorical must be TRUE or FALSE");return(NULL) }
  if(!is.na(repeatable_columns)){if(sum(repeatable_columns%in%colnames(data))!=length(repeatable_columns) ){ print("repeatable_columns must match column names");return(NULL) }}

  if(response_is_categorical==TRUE){
    family_p=switch(family_p, "poisson" = poisson(link="log"), "binomial" = binomial(link="logit"), "bionomial_log" = binomial(link="log") )
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


  nonrepeatable_columns=NULL

  for(i in 1:length(experimental_columns)){
    Data[,experimental_columns[i]]=as.factor(Data[,experimental_columns[i]])
    experimental_columns_index=c(experimental_columns_index,which(colnames(Data)==experimental_columns[i]))
    colnames(Data)[experimental_columns_index[i]]=paste("experimental_column",i,sep="")

    if(i!=1&&!experimental_columns[i]%in%repeatable_columns){
      nonrepeatable_columns=c(nonrepeatable_columns, paste("experimental_column",i,sep=""))
    }


    cat("\n")
    print(paste("_________________________________",experimental_columns[i]," is assigned to experimental_column",i,sep=""))
    cat("\n")
  }


  if(length(experimental_columns)>=2){
        for(r in 2:length(experimental_columns)){
      if(colnames(Data)[experimental_columns_index[r]]%in%nonrepeatable_columns){
        Data[,experimental_columns_index[r]]=paste(Data[,experimental_columns_index[r-1]], Data[,experimental_columns_index[r]],sep="_")
      }
    }
  }




  colnames(Data)[which(colnames(Data)==condition_column)]="condition_column"
  colnames(Data)[which(colnames(Data)==response_column)]="response_column"
  if(!is.na(covariate)) colnames(Data)[which(colnames(Data)==covariate)]="covariate"


  ####### run the formula

  if(is.na(covariate)){
    if(response_is_categorical==FALSE){
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
    if(response_is_categorical==FALSE){
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
  Data[,"condition_column"]=as.numeric(as.character(Data[,"condition_column"]))

  if(condition_is_categorical==TRUE){

    #boxplot
    boxplot( as.formula( paste0( "residual ~  condition_column") ), data=Data , xlab=condition_column, ylab=paste0(response_column, " Residual Value"), main=NULL)

    #95% CI plot
    plotData = Data %>%group_by(as.factor(condition_column))%>%dplyr::summarise(mean_residual = mean(residual), se = sd(residual)/sqrt(n()))

    colnames(plotData)[1]="condition"

    gp=ggplot2::ggplot(plotData,                ### The data frame to use.
           aes(x = condition,
               y = mean_residual)) +
      geom_errorbar(aes(ymin = mean_residual - 1.96*se,
                        ymax = mean_residual + 1.96*se),
                    width = 0.05,
                    size  = 0.5) +
      geom_point(shape = 15,
                 size  = 4) +
      theme_bw() +
      theme(axis.title   = element_text(face  = "bold")) + ylab("Mean Residual Value") + xlab("Classification")

    print(gp)

  }else{


    plot(Data[,"condition_column"], Data[,"residual"], xlab=condition_column, ylab=paste0(response_column, " Residual Value"), main=NULL)
    abline(lm( as.formula( paste0( "residual ~  condition_column") ), data=Data),col='blue')






    }

  return(Data)
}
