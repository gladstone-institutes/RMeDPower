#' @title transform_data
#'
#' @description This functions makes quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values after removing outliers. To detect outliers, the function uses Rosner's test.
#'
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of the variable related to experimental design such as "experiment", "plate", and "cell_line".
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param response_is_categorical Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param alpha numeric scalar between 0 and 1 indicating the Type I error associated with the test of outliers
#'
#' @return Quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values
#'
#' @export
#'
#' @examples transform_data(data,"classif",c("experiment","line"),"feature1","TRUE")



transform_data<-function(data,condition_column, experimental_columns, response_column, condition_is_categorical,
                         response_is_categorical=FALSE, alpha=0.05){


  data_updated=data

  ##raw qq
  check_normality(data,condition_column, experimental_columns, response_column,  condition_is_categorical,
                  response_is_categorical=FALSE, image_title="QQplot (raw data)")


  ############rosner's test begin
  trait=data[,response_column]

  options(warn=-1)
  upper_bound <- median(trait) + 3 * mad(trait)
  upper_bound
  
  
  lower_bound <- median(trait) - 3 * mad(trait)
  lower_bound


  outlierC=sum(trait>upper_bound)+sum(trait<lower_bound)
  outlierC

  if(outlierC >0  ){

    ###rosner's test
    test <- EnvStats::rosnerTest(trait,
                                 k = outlierC, alpha=alpha)

    outliers=test$all.stats$Value[test$all.stats$Outlier]

    if(length(test$all.stats$Outlier) !=0  ){
      if(! ( length(test$all.stats$Outlier) ==1 & sum(is.infinite(test$all.stats$Outlier)) ) ){

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
      }



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
                                 k = outlierC, alpha=alpha)

    if(length(test$all.stats$Outlier)>0){
      if(! ( length(test$all.stats$Outlier) ==1 & is.infinite(test$all.stats$Outlier) ) ){
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

    }

    }else{
      print("No outlier detected from the log transformed data")
    }


  }else{
    print("No outlier detected from the raw data")
  }

  return(data_updated)
}
