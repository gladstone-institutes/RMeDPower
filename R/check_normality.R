#' @title check_normality
#'
#' @description This function makes a quantile-quantile (qq) plot of the residual values of the mixed effects model. Users can check the normality of residual values by examining qqplot.
#'
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of the variable related to experimental design such as "experiment", "plate", and "cell_line".
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param repeatable_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param response_is_categorical Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param image_title title of the qq plot
#' @param na.action "complete": missing data is not allowed in all columns (default), "unique": missing data is not ollowed only in condition, experimental, and response columns
#'
#' @return A quantile-quantile (qq) plot of residual values in a mixed-effects model
#'
#' @export
#'
#' @examples check_normality(data,"classif",c("experiment","line"),"feature1","TRUE")



check_normality<-function(data, condition_column, experimental_columns, response_column,  condition_is_categorical,
                          repeatable_columns = NA, response_is_categorical=FALSE, image_title = NULL, na.action="complete"){

  if(!condition_column%in%colnames(data)){ print("condition_column should be one of the column names");return(NULL) }
  if(sum(experimental_columns%in%colnames(data))!=length(experimental_columns) ){ print("experimental_columns must match column names");return(NULL) }
  if(!response_column%in%colnames(data)){  print("response_column should be one of the column names");return(NULL) }
  if(!is.na(repeatable_columns)){if(sum(repeatable_columns%in%colnames(data))!=length(repeatable_columns) ){ print("repeatable_columns must match column names");return(NULL) }}



  if(na.action=="complete"){

    notNAindex=which( rowSums(is.na(data)) == 0 )

  }else if(na.action=="unique"){

    notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column)])) == 0 )

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
  if(response_is_categorical==TRUE) {
    cat("\n")
    print("_________________________________Categorical response variable is not accepted in the current version")}
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


  if(length(experimental_columns)>2){
        for(r in 2:length(experimental_columns)){
      if(colnames(Data)[experimental_columns_index[r]]%in%nonrepeatable_columns){
        Data[,experimental_columns_index[r]]=paste(Data[,experimental_columns_index[r-1]], Data[,experimental_columns_index[r]],sep="_")
      }
    }
  }



  colnames(Data)[which(colnames(Data)==condition_column)]="condition_column"
  colnames(Data)[which(colnames(Data)==response_column)]="response_column"



  ####### run the formula

  if(length(experimental_columns)==1){
    lmerFit <- lme4::lmer(response_column ~ condition_column + (1 | experimental_column1), data=Data)
  }else if(length(experimental_columns)==2){
    lmerFit <- lme4::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2), data=Data)
  }else if(length(experimental_columns)==3){
    lmerFit <- lme4::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3), data=Data)
  }else if(length(experimental_columns)==4){
    lmerFit <- lme4::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4), data=Data)
  }else if(length(experimental_columns)==5){
    lmerFit <- lme4::lmer(response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) + (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data)
  }

  lmerFit_s=summary(lmerFit)
  ###QQ plot

    qqnorm(lmerFit_s$residuals, ylab="Standardized Residuals", xlab="Normal Scores", main=image_title)
    qqline(lmerFit_s$residuals)

    return(lmerFit_s$residuals)



}



