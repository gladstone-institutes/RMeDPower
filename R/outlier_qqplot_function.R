#' @title outlier_qqplot_function
#'
#'
#' @description
#'
#' @param
#' @return
#'
#'
#' @export
#' @examples
#'
#library(DescTools)


paste_columns<-function(data1,data2){
  data1=as.matrix(data1)
  data2=as.matrix(data2)


   if( ncol(data1) == 1 & ncol(data2) == 1){
       if( sum(data1 == data2)==nrow(data1)){
           return(as.matrix(data1) )
       }
  }
  {
    return( as.matrix(paste(paste_columns(data1[,1:(ncol(data1)-1)],data1[,ncol(data1)]),as.matrix(data2),sep="_") ) )
  }
}

rename_experimental_variables<-function(data, experimental_columns){



  print("Experimental column names should be sorted based on a hierarchy. ex) c('Experiment', 'plate', 'cell_line') ")

  i=length(experimental_columns)

    result=NULL
    if(i==2){

      result=as.matrix( paste(data[,experimental_columns[1]], data[,experimental_columns[2]], sep="_") )

    }else{
      while(i>1){
        result=cbind(paste_columns(data[,experimental_columns[1:(i-1)]], as.matrix(data[,experimental_columns[i]]) ),result)
        i=i-1
      }

    }

  colnames(result)=paste0(experimental_columns[2:length(experimental_columns)],"_v2")
  return(cbind(result,data))


}

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




transform_data<-function(data,condition_column, experimental_columns, response_column, condition_is_categorical,
                      response_is_categorical=FALSE){


  data_updated=data

  ##raw qq
  check_normality(data,condition_column, experimental_columns, response_column,  condition_is_categorical,
                            response_is_categorical=FALSE, image_title="QQplot (raw data)")


  ############rosner's test begin
  trait=data[,response_column]

  options(warn=-1)
  upper_bound <- median(trait) + 3 * mad(trait)
  upper_bound


  outlierC=sum(trait>upper_bound)
  outlierC

  if(outlierC >0  ){

    ###rosner's test
    test <- EnvStats::rosnerTest(trait,
                       k = outlierC)

    outliers=test$all.stats$Value[test$all.stats$Outlier]

    if(length(test$all.stats$Outlier) !=0 & !is.finite(test$all.stats$Outlier)){
      cutoff=min( outliers[which(outliers>median(trait))] )


      ###plot the distribution and point the outlier boundary
      hist(trait,breaks=1000, main=paste0("Histogram of raw ",response_column, " values and detected outliers" ) )
      abline(v=cutoff,col="red")

      mtext(paste0("Cutoff at ",cutoff),cex=1.2)


      ############rosner's test end


      ########################from here, check qq after removing outliers


      data_noOutlier=data
      data_noOutlier[data_noOutlier[,response_column]>cutoff,]=NA

      data_updated=cbind(data_noOutlier[,response_column], data_updated)
      colnames(data_updated)[1]=paste0(response_column,"_noOutlier")


      check_normality(data_noOutlier,condition_column, experimental_columns, response_column,  condition_is_categorical,
                      response_is_categorical=FALSE, image_title="QQplot (outlier excluded data)")


    }else{
      print("No outlier detected from the raw data")
    }


  }else{
    print("No outlier detected from the raw data")
  }









########################do the same using log transformed values

###log transform

  data_log=data
  data_log[,response_column]=log(data_log[,response_column]+0.1^10)

  data_updated=cbind(data_log[,response_column], data_updated)
  colnames(data_updated)[1]=paste0(response_column,"_logTransformed")


  check_normality(data_log,condition_column, experimental_columns, response_column,  condition_is_categorical,
                  response_is_categorical=FALSE, image_title="QQplot (log transformed data)")



trait = data_log[,response_column] ###change the feature as you want




############rosner's test begin
options(warn=-1)
upper_bound <- median(trait) + 3 * mad(trait)
upper_bound


outlierC=sum(trait>upper_bound)
outlierC

if(outlierC >0){

###rosner's test
test <- EnvStats::rosnerTest(trait,
                   k = outlierC)

  if(length(test$all.stats$Outlier)>0 & !is.finite(test$all.stats$Outlier)){
    outliers=test$all.stats$Value[test$all.stats$Outlier]
    cutoff=min( outliers[which(outliers>median(trait))] )


    ###plot the distribution and point the outlier boundary
    hist(trait,breaks=1000, main=paste0("Histogram of log-transformed ",response_column, " values and detected outliers" ) )
    abline(v=cutoff,col="red")

    mtext(paste0("Cutoff at ",cutoff),cex=1.2)

    ############rosner's test end


    ########################from here, check qq after removing outliers


    data_log_noOutlier=data_log
    data_log_noOutlier[data_log_noOutlier[,response_column]>cutoff,]=NA
    data_updated=cbind(data_log_noOutlier[,response_column], data_updated)
    colnames(data_updated)[1]=paste0(response_column,"_logTransformed_noOutlier")

    ###qqplot
    check_normality(data_log_noOutlier,condition_column, experimental_columns, response_column,  condition_is_categorical,
                    response_is_categorical=FALSE, image_title="QQplot (log transformed & ouliter excluded data)")




  }else{
    print("No outlier detected from the log transformed data")
  }


}else{
  print("No outlier detected from the raw data")
}

return(data_updated)
}
