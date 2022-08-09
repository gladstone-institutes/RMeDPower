#' @title transform_data_by_residual
#'
#' @description This functions makes quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values after removing outliers. To detect outliers, the function uses Rosner's test.
#'
#'
#' @param data Input data
#' @param residual_values User provided residual values
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param alpha numeric scalar between 0 and 1 indicating the Type I error associated with the test of outliers
#'
#' @return Quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values
#'
#' @export
#'
#' @examples transform_data(data,"classif",c("experiment","line"),"feature1","TRUE")



transform_data_by_residual <-function(data, residual_values, response_column,  alpha=0.05){



  if(length(residual_values)!=nrow(data)){
    print("Error: The two input data have different number of rows.")
  }
  if(!response_column%in%colnames(data)){  print("response_column should be one of the column names");return(NULL) }

  Data_updated=cbind(residual_values, data)
  colnames(Data_updated)[1]="residuals"

  Data_updated=as.data.frame(Data_updated)

  ##raw qq
  qqnorm(Data_updated$residuals, ylab="Standardized Residuals", xlab="Normal Scores", main="QQplot (raw data)")
  qqline(Data_updated$residuals)



  ############rosner's test begin
  trait=Data_updated$residuals

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
          Data_noOutlier[Data_updated$residuals>cutoff1,]=NA
        }
        if(!is.na(cutoff2)){
          Data_noOutlier[Data_updated$residuals<cutoff2,]=NA
        }


        Data_updated=cbind(Data_noOutlier[,response_column], Data_updated)
        colnames(Data_updated)[1]=paste0(response_column,"_noOutlier")


        ##qqplot
        qqnorm(Data_updated$residuals, ylab="Standardized Residuals", xlab="Normal Scores", main="QQplot (outlier excluded Data))")
        qqline(Data_updated$residuals)

      }



    }else{
      print("No outlier detected from the raw Data")
    }


  }else{
    print("No outlier detected from the raw Data")
  }








  return(Data_updated)
}
