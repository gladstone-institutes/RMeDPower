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
#' @param repeatable_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param response_is_categorical Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param alpha numeric scalar between 0 and 1 indicating the Type I error associated with the test of outliers
#' @param na.action "complete": missing data is not allowed in all columns (default), "unique": missing data is not ollowed only in condition, experimental, and response columns
#'
#' @return Quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values
#'
#' @export
#'
#' @examples transform_data(data,"classif",c("experiment","line"),"feature1","TRUE")



transform_data<-function(data, condition_column, experimental_columns, response_column, condition_is_categorical,
                         repeatable_columns = NA, response_is_categorical=FALSE, alpha=0.05, na.action="complete"){


  if(na.action=="complete"){

    notNAindex=which( rowSums(is.na(data)) == 0 )

  }else if(na.action=="unique"){

    notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column)])) == 0 )

  }



  data=data[notNAindex,]


  Data_updated=data

  ##raw qq
  residual=check_normality(data, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                  response_column = response_column,  condition_is_categorical = condition_is_categorical,
                  response_is_categorical = FALSE, image_title="QQplot (raw data)", na.action=na.action)




  ############rosner's test begin
  trait=residual

  options(warn=-1)
  upper_bound <- median(trait) + 3 * mad(trait)
  upper_bound

  lower_bound <- median(trait) - 3 * mad(trait)
  lower_bound

  outlierC=sum(trait>upper_bound)+sum(trait<lower_bound)
  outlierC

  if( outlierC >0 ){

    ###rosner's test
    test <- EnvStats::rosnerTest(trait,
                                 k = outlierC, alpha=alpha)

    outliers=test$all.stats$Value[test$all.stats$Outlier]

    if(length(test$all.stats$Outlier) !=0  ){
      if(! ( length(test$all.stats$Outlier) ==1 &  sum(is.infinite(test$all.stats$Outlier)) ) ) {

        cutoff1=NA
        cutoff2=NA
        if(sum(outliers>median(trait)) >0) cutoff1=min( outliers[which(outliers>median(trait))] )
        if(sum(outliers<median(trait)) >0) cutoff2=max( outliers[which(outliers<median(trait))] )


        ###plot the distribution and point the outlier boundary
        hist(trait,breaks=1000, main=paste0("Histogram of raw ",response_column, " values and detected outliers" ) )

        if(!is.na(cutoff1)) abline(v=cutoff1,col="red")
        if(!is.na(cutoff2)) abline(v=cutoff2,col="red")

        if(!is.na(cutoff1)&is.na(cutoff2)){
          mtext(paste0("Cutoff at ",cutoff1 ),cex=1.2)
        }else if(!is.na(cutoff2)&is.na(cutoff1)){
          mtext(paste0("Cutoff at ",cutoff2 ),cex=1.2)
        }else{
          mtext(paste0("Cutoff at ",cutoff1, " and ", cutoff2 ),cex=1.2)
        }



        ############rosner's test end


        ########################from here, check qq after removing outliers


        Data_noOutlier=Data_updated
        if(!is.na(cutoff1)){
          Data_noOutlier[residual>cutoff1,]=NA
        }
        if(!is.na(cutoff2)){
          Data_noOutlier[residual<cutoff2,]=NA
        }


        Data_updated=cbind(Data_noOutlier[,response_column], Data_updated)
        colnames(Data_updated)[1]=paste0(response_column,"_noOutlier")


        check_normality(Data_noOutlier, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                        response_column = response_column,  condition_is_categorical = condition_is_categorical,
                        response_is_categorical = FALSE,  image_title="QQplot (outlier excluded Data)", na.action=na.action)
      }



    }else{
      print("No outlier detected from the raw Data")
    }


  }else{
    print("No outlier detected from the raw Data")
  }









  ########################do the same using log transformed values

  ###log transform

  Data_log=data
  temp1=Data_log[,response_column]
  temp2=min(temp1[temp1>0], temp1[temp1<0])
  temp2 = 1 - min(temp1)
  Data_log[,response_column]=log(Data_log[,response_column]+temp2)

  Data_updated=cbind(Data_log[,response_column], Data_updated)
  colnames(Data_updated)[1]=paste0(response_column,"_logTransformed")


  residual=check_normality(Data_log, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                  response_column = response_column,  condition_is_categorical = condition_is_categorical,
                  response_is_categorical = FALSE,  image_title="QQplot (log transformed Data)", na.action=na.action)



  trait = residual ###change the feature as you want




  ############rosner's test begin
  options(warn=-1)
  upper_bound <- median(trait) + 3 * mad(trait)
  upper_bound

  lower_bound <- median(trait) - 3 * mad(trait)
  lower_bound

  outlierC=sum(trait>upper_bound)+sum(trait<lower_bound)
  outlierC

  if(outlierC >0){

    ###rosner's test
    test <- EnvStats::rosnerTest(trait,
                                 k = outlierC, alpha=alpha)

    if(length(test$all.stats$Outlier)>0){
      if(! ( length(test$all.stats$Outlier) ==1 &  sum(is.infinite(test$all.stats$Outlier)) ) ) {
        outliers=test$all.stats$Value[test$all.stats$Outlier]
        cutoff1=NA
        cutoff2=NA
        if(sum(outliers>median(trait)) >0) cutoff1=min( outliers[which(outliers>median(trait))] )
        if(sum(outliers<median(trait)) >0) cutoff2=max( outliers[which(outliers<median(trait))] )


        ###plot the distribution and point the outlier boundary
        hist(trait,breaks=1000, main=paste0("Histogram of log-transformed ",response_column, " values and detected outliers" ) )
        if(!is.na(cutoff1)) abline(v=cutoff1,col="red")
        if(!is.na(cutoff2)) abline(v=cutoff2,col="red")

        if(!is.na(cutoff1)&is.na(cutoff2)){
          mtext(paste0("Cutoff at ",cutoff1 ),cex=1.2)
        }else if(!is.na(cutoff2)&is.na(cutoff1)){
          mtext(paste0("Cutoff at ",cutoff2 ),cex=1.2)
        }else{
          mtext(paste0("Cutoff at ",cutoff1, " and ", cutoff2 ),cex=1.2)
        }


        ############rosner's test end


        ########################from here, check qq after removing outliers


        Data_log_noOutlier=Data_log

        if(!is.na(cutoff1)){
          Data_log_noOutlier[residual>cutoff1,]=NA
        }
        if(!is.na(cutoff2)){
          Data_log_noOutlier[residual<cutoff2,]=NA
        }

        Data_updated=cbind(Data_log_noOutlier[,response_column], Data_updated)
        colnames(Data_updated)[1]=paste0(response_column,"_logTransformed_noOutlier")

        ###qqplot
        check_normality(Data_log_noOutlier, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                        response_column = response_column,  condition_is_categorical = condition_is_categorical,
                        response_is_categorical = FALSE,  image_title="QQplot (log transformed & ouliter excluded Data)", na.action=na.action)

    }

    }else{
      print("No outlier detected from the log transformed Data")
    }


  }else{
    print("No outlier detected from the raw Data")
  }


  return(Data_updated)
}
