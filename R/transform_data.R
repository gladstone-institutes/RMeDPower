#' @title transform_data
#'
#' @export


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
