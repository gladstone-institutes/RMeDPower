#' @title check_normality
#'
#' @export



check_normality<-function(data,condition_column, experimental_columns, response_column,  condition_is_categorical,
                          response_is_categorical=FALSE, image_title=""){

  if(!condition_column%in%colnames(data)){ print("condition_column should be one of the column names");return(NULL) }
  if(sum(experimental_columns%in%colnames(data))!=length(experimental_columns) ){ print("experimental_columns must match column names");return(NULL) }
  if(!response_column%in%colnames(data)){  print("response_column should be one of the column names");return(NULL) }


  ####### remove empty lines
  Data <- data[complete.cases(data),]

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


  for(i in 1:length(experimental_columns)){
    Data[,experimental_columns[i]]=as.factor(Data[,experimental_columns[i]])
    experimental_columns_index=c(experimental_columns_index,which(colnames(Data)==experimental_columns[i]))
    colnames(Data)[experimental_columns_index[i]]=paste("experimental_column",i,sep="")
    cat("\n")
    print(paste("_________________________________",experimental_columns[i]," is assigned to experimental_column",i,sep=""))
    cat("\n")
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

}



