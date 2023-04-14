#' @title transform_data2_covariate
#'
#' @description This functions makes quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values after removing outliers. To detect outliers, the function uses Rosner's test.
#'
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of the variable related to experimental design such as "experiment", "plate", and "cell_line". They should be in order, for example, "experiment" should always come first .
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param covariate The name of the covariate to control in the regression model
#' @param method The method used to detect outliers. "rosner" (default) runs Rosner's test and "cook" runs Cook's distance.
#' @param repeatable_columns Name of experimental variables that may appear repeatedly with the same ID. For example, cell_line C1 may appear in multiple experiments, but plate P1 cannot appear in more than one experiment
#' @param response_is_categorical Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param family_p The type of distribution family to specify when the response is categorical. If family is "binary" then binary(link="log") is used, if family is "poisson" then poisson(link="logit") is used, if family is "poisson_log" then poisson(link=") log") is used.
#' @param alpha numeric scalar between 0 and 1 indicating the Type I error associated with the test of outliers
#' @param na.action "complete": missing data is not allowed in all columns (default), "unique": missing data is not allowed only in condition, experimental, and response columns. Selecting "complete" removes an entire row when there is one or more missing values, which may affect the distribution of other features.
#'
#' @return For continuous data, the function returns quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values. For discrete data, it returns a histogram. If "rosner" is chosen, a matrix with updated feature values after transformation will be returned. If "cook" is choose, a list with  a matrix with updated feature values after transformation will be returned,
#'
#' @export
#'
#' @examples result=transform_data2(data=data, condition_column="classif", experimental_columns=c("experiment","line"), response_column="feature", condition_is_categorical=TRUE, response_is_categorical=FALSE, alpha=0.05, repeatable_columns = "line", method="cook", na.action="complete")
#' @examples result=transform_data2(data=data, condition_column="classif", experimental_columns=c("experiment","line"), response_column="feature", condition_is_categorical=TRUE, response_is_categorical=FALSE, alpha=0.05, repeatable_columns = "line", method="cook", na.action="complete")
#' @examples result=transform_data2(data=data, condition_column="classif", experimental_columns=c("experiment","line"), response_column="feature", condition_is_categorical=TRUE, response_is_categorical=TRUE, family_p="poisson", alpha=0.05, repeatable_columns = "line", method="cook", na.action="complete")


rosner_test<- function (trait, response_column, alpha, hist_text) {

  cutoff1=NA
  cutoff2=NA

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


        if(sum(outliers>median(trait)) >0) cutoff1=min( outliers[which(outliers>median(trait))] )
        if(sum(outliers<median(trait)) >0) cutoff2=max( outliers[which(outliers<median(trait))] )


        ###plot the distribution and point the outlier boundary
        hist(trait,breaks=1000, main=paste0("Histogram of ", hist_text, " ", response_column, " values and detected outliers" ) )

        if(!is.na(cutoff1)) abline(v=cutoff1,col="red")
        if(!is.na(cutoff2)) abline(v=cutoff2,col="red")

        if(!is.na(cutoff1)&is.na(cutoff2)){
          mtext(paste0("Cutoff at ",cutoff1 ),cex=1.2)
        }else if(!is.na(cutoff2)&is.na(cutoff1)){
          mtext(paste0("Cutoff at ",cutoff2 ),cex=1.2)
        }else if(!is.na(cutoff1)&!is.na(cutoff2)){
          mtext(paste0("Cutoff at ",cutoff2, " and ", cutoff1 ),cex=1.2)
        }



        ############rosner's test end


      }



    }else{
      print("No outlier detected from the raw Data")
    }


  }else{
    print("No outlier detected from the raw Data")
  }

  return(c(cutoff2,cutoff1))
}

 # model=lms[[1]]
 # Data=lms[[2]]
 # experimental_columns=lms[[3]]
 # residuals=lms[[4]]
 # response_column=response_column
 # hist_text="raw"
 # response_is_categorical=response_is_categorical

cooks_test<- function (model, fixed_global_variable_data, experimental_columns, residuals, response_column, hist_text) {


  cooks_result=lapply(1:length(experimental_columns),
         function(i){

             alt.est <- influence.ME::influence(model, group=experimental_columns[i])


           cooks.distance(alt.est)
         }
  )


  alt.est <- influence.ME::influence(model, obs=TRUE)


  cooks_result_sample=cooks.distance(alt.est)
  cooks_result=c(list(cooks_result_sample), cooks_result)
  names(cooks_result)=c("cooks_distance_sample", paste0("cooks_distance_",experimental_columns) )

  filtering_targets=cooks_result_sample>4/nrow(fixed_global_variable_data)


  conditions=unique(fixed_global_variable_data$condition_column)

  lapply(conditions, function(conditions_p)
  {
      outliers=fixed_global_variable_data[filtering_targets & fixed_global_variable_data$condition_column%in%conditions_p, "response_column"]
      Data_temp=fixed_global_variable_data[fixed_global_variable_data$condition_column%in%conditions_p, "response_column"]
      #residual_temp=residuals[Data$condition_column%in%conditions_p]

      #plot histogram

      cutoff1=NA
      cutoff2=NA

      medianV=median(Data_temp)
      if(sum(outliers>medianV) >0) cutoff1=min( outliers[which(outliers>medianV)] )
      if(sum(outliers<medianV) >0) cutoff2=max( outliers[which(outliers<medianV)] )


      ###plot the distribution and point the outlier boundary
      hist(Data_temp, breaks=1000,
           main=paste0("Histogram of ", hist_text," ", response_column, " values and detected outliers in condition ",conditions_p ) )

      if(!is.na(cutoff1)) abline(v=cutoff1,col="red")
      if(!is.na(cutoff2)) abline(v=cutoff2,col="red")

      if(!is.na(cutoff1)&is.na(cutoff2)){
        mtext(paste0("Cutoff at ",cutoff1 ),cex=1.2)
      }else if(!is.na(cutoff2)&is.na(cutoff1)){
        mtext(paste0("Cutoff at ",cutoff2 ),cex=1.2)
      }else if(!is.na(cutoff1)&!is.na(cutoff2)){
        mtext(paste0("Cutoff at ",cutoff2, " and ", cutoff1 ),cex=1.2)
      }
    }
  )



  return(cooks_result)

}


transform_data2_covariate<-function(data, condition_column, experimental_columns, response_column, condition_is_categorical, covariate=NA, method="rosner",
                         repeatable_columns = NA, response_is_categorical=FALSE, family_p=NULL, alpha=0.05, na.action="complete"){


  if(na.action=="complete"){

    notNAindex=which( rowSums(is.na(data)) == 0 )

  }else if(na.action=="unique"){

    if(is.na(covariate)) notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column, covariate)])) == 0 )
    else notNAindex=which( rowSums(is.na(data[,c(condition_column, experimental_columns, response_column)])) == 0 )


  }



  data=data[notNAindex,]


  Data_updated=data

  if(response_is_categorical==FALSE){
    ##raw qq
    residual=check_normality_covariate(data, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                             response_column = response_column,  condition_is_categorical = condition_is_categorical, covariate=covariate,
                             response_is_categorical = response_is_categorical, image_title="QQplot (raw data)", na.action=na.action)

  }




  if(method=="rosner"){

    cutoffs=rosner_test(trait=residual, response_column=response_column, alpha=alpha, hist_text = "raw residual")

    if(sum(is.na(cutoffs))<2){

      Data_noOutlier=Data_updated
      if(!is.na(cutoffs[1])){
        Data_noOutlier[residual<=cutoffs[1],]=NA
      }
      if(!is.na(cutoffs[2])){
        Data_noOutlier[residual>=cutoffs[2],]=NA
      }


      Data_updated=cbind(Data_noOutlier[,response_column], Data_updated)
      colnames(Data_updated)[1]=paste0(response_column,"_noOutlier")


      check_normality_covariate(Data_noOutlier, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                      response_column = response_column,  condition_is_categorical = condition_is_categorical, covariate=covariate,
                      response_is_categorical = FALSE,  image_title="QQplot (outlier excluded Data)", na.action=na.action)


    }

  }else if(method=="cook"){


    #run regression
    lms=get_model_and_data_covariate(data=data, condition_column=condition_column, experimental_columns=experimental_columns,
                                 response_column=response_column, condition_is_categorical=condition_is_categorical, covariate=covariate,
                                 repeatable_columns=repeatable_columns, response_is_categorical=response_is_categorical, family_p=family_p, na.action=na.action)

    #run cook
    fixed_global_variable_data<<-lms[[2]]
    family_p<<-family_p

    cooks_result=cooks_test(lms[[1]], lms[[2]], lms[[3]], lms[[4]], response_column=response_column, hist_text="raw")

    rm(fixed_global_variable_data)
    rm(family_p)
    #get cutoff

    filtering_targets=cooks_result[[1]]>4/nrow(data)

    if(sum(filtering_targets)>0){

      Data_noOutlier=Data_updated
      Data_noOutlier[filtering_targets,]=NA



      Data_updated=cbind(Data_noOutlier[,response_column], Data_updated)
      colnames(Data_updated)[1]=paste0(response_column,"_noOutlier")

      if(response_is_categorical==FALSE){
        check_normality_covariate(Data_noOutlier, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                        response_column = response_column,  condition_is_categorical = condition_is_categorical, covariate=covariate,
                        response_is_categorical = FALSE,  image_title="QQplot (outlier excluded Data)", na.action=na.action)

      }



    }


  }









  if(response_is_categorical==FALSE){

    ########################do the same using log transformed values

    ###log transform

    Data_log=data
    temp1=Data_log[, response_column]


    #seperate zero and negative
    if(min(Data_log[,response_column])<0){

      temp2 = abs(min(temp1))/10 - min(temp1)
      Data_log[,response_column]=log(Data_log[,response_column]+temp2)

    }else if(min(Data_log[,response_column])==0){

      temp2 = min(temp1[temp1>0])/10
      Data_log[,response_column]=log(Data_log[,response_column]+temp2)

    }else{

      Data_log[,response_column]=log(Data_log[,response_column])

    }


    Data_updated=cbind(Data_log[,response_column], Data_updated)
    colnames(Data_updated)[1]=paste0(response_column,"_logTransformed")


    residual=check_normality_covariate(Data_log, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                             response_column = response_column,  condition_is_categorical = condition_is_categorical, covariate=covariate,
                             response_is_categorical = FALSE,  image_title="QQplot (log transformed Data)", na.action=na.action)



    trait = residual ###change the feature as you want



    if(method=="rosner"){

      cutoffs=rosner_test(trait=residual, response_column=response_column, alpha=alpha, hist_text = "log residual")

      if(sum(is.na(cutoffs))<2){


        Data_log_noOutlier=Data_log

        if(!is.na(cutoffs[1])){
          Data_log_noOutlier[residual<=cutoffs[1],]=NA
        }
        if(!is.na(cutoffs[2])){
          Data_log_noOutlier[residual>=cutoffs[2],]=NA
        }

        Data_updated=cbind(Data_log_noOutlier[,response_column], Data_updated)
        colnames(Data_updated)[1]=paste0(response_column,"_logTransformed_noOutlier")

        ###qqplot
        check_normality_covariate(Data_log_noOutlier, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                        response_column = response_column,  condition_is_categorical = condition_is_categorical, covariate=covariate,
                        response_is_categorical = FALSE,  image_title="QQplot (log transformed & ouliter excluded Data)", na.action=na.action)

      }

    }else if(method=="cook"){

      #run regression
      lms=get_model_and_data_covariate(data=Data_log, condition_column=condition_column, experimental_columns=experimental_columns,
                             response_column=response_column, condition_is_categorical=condition_is_categorical, covariate=covariate,
                             repeatable_columns=repeatable_columns, response_is_categorical=response_is_categorical, family_p=family_p, na.action=na.action)

      #run cook
      fixed_global_variable_data<<-lms[[2]]
      family_p<<-family_p

      cooks_result2=cooks_test(lms[[1]],lms[[2]], lms[[3]], lms[[4]], response_column=response_column, hist_text="log transformed")

      #get cutoff

      rm(fixed_global_variable_data)
      rm(family_p)

      filtering_targets=cooks_result2[[1]]>4/nrow(data)

      if(sum(filtering_targets)>0){

        Data_log_noOutlier=Data_log
        Data_log_noOutlier[filtering_targets,]=NA


        Data_updated=cbind(Data_log_noOutlier[,response_column], Data_updated)
        colnames(Data_updated)[1]=paste0(response_column,"_logTransformed_noOutlier")


        check_normality_covariate(Data_log_noOutlier, condition_column = condition_column, experimental_columns = experimental_columns,  repeatable_columns = repeatable_columns,
                        response_column = response_column,  condition_is_categorical = condition_is_categorical, covariate=covariate,
                        response_is_categorical = FALSE,  image_title="QQplot (log transformed & outlier excluded Data)", na.action=na.action)


      }

    }






    if(method=="rosner"){

      return(Data_updated)

    }else if(method=="cook"){
      result=c(list(Data_updated ), cooks_result, cooks_result2)
      names(result)[1]="Data_updated"
      names(result)[5:7]=paste0(names(result)[5:7],"_logTransformed")

      return(result)
    }

  }else{

    if(method=="rosner"){

      return(Data_updated)

    }else if(method=="cook"){
      result=c(list(Data_updated ), cooks_result)
      names(result)[1]="Data_updated"
      return(result)
    }

  }




}
