#' @title calculate_power
#'
#'
#' @description This function uses simulation to perform power analysis. It is designed to explore the power of biological experiments and to suggest an optimal number of experimental variables
#' @description with reasonable power. The backbone of the function is based on simr package, which fits a fixed effect or mixed effect model based on the observed data and simulates response
#' @description variables. Users can test the power of different combinations of experimental variables and parameters.
#'
#' Note: The current version does not accept categorical response variables, sample size parameters smaller than the observed samples size
#'
#' @import multtest
#' @import simr
#' @import lme4
#' @import lmerTest
#' @import readxl
#'
#' @param data Input data
#' @param condition_column Name of the condition variable (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_columns Name of the variable related to experimental design such as "experiment", "plate", and "cell_line".
#' @param response_column Name of the variable observed by performing the experiment. ex) intensity.
#' @param power_curve 1: Power simulation over a range of sample sizes or levels. 0: Power calculation over a single sample size or a level.
#' @param condition_is_categorical Specify whether the condition variable is categorical. TRUE: Categorical, FALSE: Continuous.
#' @param response_is_categorical Default: the observed variable is continuous Categorical response variable will be implemented in the future. TRUE: Categorical , FALSE: Continuous (default).
#' @param nsim The number of simulations to run. Default=1000
#' @param target_columns Name of the experimental parameters to use for the power calculation.
#' @param levels 1: Amplify the number of corresponding target parameter. 0: Amplify the number of samples from the corresponding target parameter, ex) If target_columns = c("experiment","cell_line") and if you want to expand the number of experiment and sample more cells from each cell line, levels = c(1,0).
#' @param max_size Maximum levels or sample sizes to test. Default: the current level or the current sample size x 5. ex) If max_levels = c(10,5), it will test upto 10 experiments and 5 cell lines.
#' @param breaks Levels /sample sizes of the variable to be specified along the power curve.. Default: max(1, round( the number of current levels / 5 ))
#' @param output Output file name
##### If variance estimates should be estimated from data
#' @param  effect_size If you know the effect size of your condition variable, provided it. If the effect size is not provided, it will be estimated from your data
##### If variance estimates are to be assigned by a user
#' @param  ICC Intra-Class Coefficients (ICC) for each parameter
#' @return A power curve image or a power calculation result printed in a text file
#'
#' @export
#' @examples

calculate_power <- function(data, condition_column, experimental_columns, response_column, target_columns, power_curve, condition_is_categorical,
                            response_is_categorical=FALSE, nsimn=1000,
                            levels=NULL, max_size=NULL, breaks=NULL, effect_size=NULL, ICC=NULL, output=NULL){






  ######input error handler
  if(length(levels)!=length(target_columns)){ print("User should specify levels of all target parameters") }
  if(length(power_curve)==0 | !power_curve%in%c(0,1)){ print("power_curve must be 0 or 1");return(NULL) }
  if(!condition_column%in%colnames(data)){ print("condition_column should be one of the column names");return(NULL) }
  if(sum(experimental_columns%in%colnames(data))!=length(experimental_columns) ){ print("experimental_columns must match column names");return(NULL) }
  if(!response_column%in%colnames(data)){  print("response_column should be one of the column names");return(NULL) }
  if(response_is_categorical==TRUE){ print("response_is_categorical is TRUE. Categorical response variable is not accepted in the current version");return(NULL) }
  if(is.null(condition_is_categorical) | !condition_is_categorical%in%c(TRUE,FALSE)){ print("condition_is_categorical must be TRUE or FALSE");return(NULL) }
  if(! (is.numeric(nsimn)&&nsimn>0) ){ print("nsimn should be a positive integer");return(NULL) }
  if(sum(target_columns%in%colnames(data))!=length(target_columns) ){ print("target_columns must match column names");return(NULL) }
  if(!levels%in%c(0,1)){ print("levels must be 0 or 1");return(NULL) }
  if(!( is.null(max_size) | (is.numeric(max_size)&&sum(max_size>0)==length(max_size)) ) ){print("max_size a positive integer");return(NULL) }
  if(!( is.null(breaks) | (is.numeric(breaks)&&breaks>0) ) ){ print("breaks must be a positive integer");return(NULL) }
  if(!( is.null(effect_size) | (is.numeric(effect_size)&&effect_size>0) ) ){ print("effect_size a positive integer");return(NULL) }




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




  ###### indices of target parameters in experimental variables
  target_i=NULL
  target_columns_renamed=NULL
  ###### match target parameters
  for(i in 1:length(target_columns)){
    cn=colnames(Data)[which(colnames_original==target_columns[i])]
    target_columns_renamed[i]=cn

    target_i=c(target_i,as.integer(substr(cn,nchar(cn),nchar(cn))))
  }




  ####### run the formula
  if(length(ICC)==0){
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
  }else{

    lmerFit=stats::lm(response_column ~ condition_column, data=Data)

  }





  slmerFit <- summary(lmerFit)
  cat("\n")
  print("__________________________________________________________________Model statistics:")
  print(slmerFit)
  cat("\n")

  fixed_effects=slmerFit$coefficients[,1]





  ##### ICC based variance estimation
  if(length(ICC)>0){

    ####### If there is only one category, add one more
    for(i in 1:length(experimental_columns)){
      if(length(table(Data[,experimental_columns_index[i]]))==1){
        categroy_temp=paste0(Data[,experimental_columns_index[i]][1],"_2")
        levels(Data[,experimental_columns_index[i]])=c(levels(Data[,experimental_columns_index[i]]),categroy_temp)
        Data[,experimental_columns_index[i]][sample(1:length(Data[,1]),round(length(Data[,1])/2))]=rep(categroy_temp,round(length(Data[,1])/2))
      }
    }




    ####### Estimated variance of the experimental variable
    if(length(experimental_columns)==1){
      varEs=(slmerFit$sigma)^2*ICC/(1-ICC)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variable:",varEs))
      cat("\n")
    }else if(length(experimental_columns)==2){

      a <- matrix(c( 1-ICC[1], -ICC[2],-ICC[1],1-ICC[2]), nrow=2, ncol=2)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2 ) , nrow=2, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      artificial_glmer=simr::makeLmer(formula = response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2)
                                      , data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )

    }else if(length(experimental_columns)==3){

      a <- matrix(c( 1-ICC[1], -ICC[2], -ICC[3], -ICC[1], 1-ICC[2], -ICC[3], -ICC[1], -ICC[2], 1-ICC[3]), nrow=3, ncol=3)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2, ICC[3]*(slmerFit$sigma)^2 ) , nrow=3, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")


      artificial_glmer=simr::makeLmer(formula = response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) +
                                        (1 | experimental_column3) , data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )


    }else if(length(experimental_columns)==4){

      a <- matrix(c( 1-ICC[1], -ICC[2], -ICC[3], -ICC[4], -ICC[1], 1-ICC[2], -ICC[3], -ICC[4], -ICC[1], -ICC[2], 1-ICC[3], -ICC[4], -ICC[1], -ICC[2], -ICC[3], 1-ICC[4]), nrow=4, ncol=4)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2, ICC[3]*(slmerFit$sigma)^2, ICC[4]*(slmerFit$sigma)^2 ) , nrow=4, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      artificial_glmer=simr::makeLmer(formula = response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) +
                                        (1 | experimental_column3) + (1 | experimental_column4) , data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )

    }else if(length(experimental_columns)==5){

      a <- matrix(c( 1-ICC[1], -ICC[2], -ICC[3], -ICC[4], -ICC[5], -ICC[1], 1-ICC[2], -ICC[3], -ICC[4], -ICC[5], -ICC[1], -ICC[2], 1-ICC[3], -ICC[4], -ICC[5]
                     , -ICC[1], -ICC[2], -ICC[3], 1-ICC[4], -ICC[5], -ICC[1], -ICC[2], -ICC[3], -ICC[4], 1-ICC[5]), nrow=4, ncol=4)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2, ICC[3]*(slmerFit$sigma)^2, ICC[4]*(slmerFit$sigma)^2, ICC[5]*(slmerFit$sigma)^2 ) , nrow=4, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      artificial_glmer=simr::makeLmer(formula = response_column ~ condition_column + (1 | experimental_column1) + (1 | experimental_column2) +
                                        (1 | experimental_column3) + (1 | experimental_column4) + (1 | experimental_column5), data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )

    }


    cat("\n")
    print(artificial_glmer)
    print("__________________________________________________________________Model statistics:")
    print(summary(artificial_glmer))
    cat("\n")
    lmerFit=artificial_glmer


  }




  ####### print observed levels and sample sizes
  maxs=NULL
  mins=NULL
  lens=NULL
  for(i in 1:length(experimental_columns)){
    cat("\n")
    print(paste("__________________________________________________________________Levels and sample sizes of",experimental_columns[i]))
    cat("\n")
    xtabs_s=stats::xtabs(~Data[,paste0("experimental_column",i)])


    print(xtabs_s)
    maxs=c(maxs,max(xtabs_s))
    if(xtabs_s[1]==0){
      mins=c(mins,min(xtabs_s[-1]))
    }else{
      mins=c(mins,min(xtabs_s))
    }

    lens=c(lens,length(xtabs_s))
    cat("\n")
    print(paste("_________________________________Max sample size:",maxs[i]))
    print(paste("_________________________________Min sample size:",mins[i]))
    print(paste("_________________________________Count of levels:",lens[i]))
    cat("\n")
  }







  ##### Assign known effect sizes
  if(length(effect_size)>0){


    fixef(lmerFit)["condition_column1"] <- effect_size

    cat("\n")
    print(paste0("_________________________________Effect size of the condition_column is now ",effect_size))
    cat("\n")

  }




  if(power_curve==0){




    ####### extend parameter levels and sample sizes


    if(length(max_size)==0){max_size=rep(0, length(target_columns))}
    if(length(breaks)==0){
      breaks=rep(0,length(target_columns))
    }


    for(i in 1:length(target_columns)){#target_columns_renamed, levels and max_size follow user input order but lens, mins, maxs don't
      if(levels[i]==1){ # increase levels

        if(max_size[i]==0){
          max_size[i]=lens[target_i[i]]*5
        }else if(max_size[i]<lens[target_i[i]]){
          cat("\n")
          print(paste("_________________________________Max size is set to ",max_size[i]," which is smaller than the observed max size ",lens[target_i[i]],". The observed max size will be used instead.",sep=""))
          cat("\n")
          max_size[i]=lens[target_i[i]]
        }
        if(breaks[i]==0) breaks[i] = max(1, round( lens[target_i[i]] / 5 ))
        if(i==1){
          extended_target_columns=simr::extend(lmerFit,along=target_columns_renamed[i],n=max_size[i])
        }else{
          extended_target_columns=simr::extend(extended_target_columns,along=target_columns_renamed[i],n=max_size[i])
        }

      }else{ # increase sample sizes
        if(max_size[i]==0){
          max_size[i]=maxs[target_i[i]]*5
        }else if(max_size[i]<maxs[target_i[i]]){
          cat("\n")
          print(paste("_________________________________Max size is set to ",max_size[i]," which is smaller than the observed max size ",maxs[target_i[i]],". The observed max size will be used instead.",sep=""))
          cat("\n")
          max_size[i]=maxs[target_i[i]]
        }
        if(breaks[i]==0) breaks[i] = max(1, round( maxs[target_i[i]] / 5 ))
        if(i==1){
          extended_target_columns=simr::extend(lmerFit,within=target_columns_renamed[i],n=max_size[i])
        }else{
          extended_target_columns=simr::extend(extended_target_columns,within=target_columns_renamed[i],n=max_size[i])
        }
      }
      cat("\n")
      print(paste("Input max size of",target_columns[i]))
      print(max_size[i])
      cat("\n")

    }

    cat("\n")
    print("__________________________________________________________________Extended parameters:")
    for(i in 1:length(target_columns)){
      if(target_columns_renamed[i]=="experimental_column1"){
        print(xtabs(~experimental_column1,data=attributes(extended_target_columns)$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column2"){
        print(xtabs(~experimental_column2,data=attributes(extended_target_columns)$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column3"){
        print(xtabs(~experimental_column3,data=attributes(extended_target_columns)$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column4"){
        print(xtabs(~experimental_column4,data=attributes(extended_target_columns)$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column5"){
        print(xtabs(~experimental_column5,data=attributes(extended_target_columns)$newData))
        cat("\n")
      }

    }
    ###### power simulation

    ps=simr::powerSim(extended_target_columns, test=simr::fixed("condition_column"),nsim=nsimn)
    cat("\n")
    print("__________________________________________________________________Power simulation result:")
    print(ps)
    cat("\n")
    if(length(output)!=0){
      sink(output)
    }else{
      sink(paste("Power_curve_experiment_",target_columns[i],"_",Sys.Date(),".txt",sep=""))
    }

    print(ps)
    sink()





  }else{#power curve =1


    ####### extend parameter levels and sample sizes
    extended_target_columns=list()
    if(length(max_size)==0){max_size=rep(0,length(target_columns))}
    if(length(breaks)==0){
      breaks=rep(0,length(target_columns))
    }





    for(i in 1:length(target_columns)){
      if(levels[i]==1){ # increase levels
        if(max_size[i]==0){
          max_size[i]=lens[target_i[i]]*5
        }else if(max_size[i]<lens[target_i[i]]){
          cat("\n")
          print(paste("_________________________________Max size is set to ",max_size[i]," which is smaller than the observed max size ",lens[target_i[i]],". The observed max size will be used instead.",sep=""))
          cat("\n")
          max_size[i]=lens[target_i[i]]
        }
        if(breaks[i]==0) breaks[i] = max(1, round( lens[target_i[i]] / 5 ))
        extended_target_columns=c(extended_target_columns,list(simr::extend(lmerFit,along=target_columns_renamed[i],n=max_size[i])))
      }else{ # increase sample sizes
        if(max_size[i]==0){
          max_size[i]=maxs[target_i[i]]*5
        }
        if(breaks[i]==0) breaks[i] = max(1, round( maxs[target_i[i]] / 5 ))
        extended_target_columns=c(extended_target_columns,list(simr::extend(lmerFit,within=target_columns_renamed[i],n=max_size[i])))
      }
      cat("\n")
      print(paste("_________________________________Power simulation will be performed based on the max size of",target_columns[i],":"))
      print(max_size[i])
      cat("\n")

    }


    cat("\n")
    print("__________________________________________________________________Extended parameters:")
    for(i in 1:length(target_columns)){
      if(target_columns_renamed[i]=="experimental_column1"){
        print(xtabs(~experimental_column1,data=attributes(extended_target_columns[[i]])$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column2"){
        print(xtabs(~experimental_column2,data=attributes(extended_target_columns[[i]])$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column3"){
        print(xtabs(~experimental_column3,data=attributes(extended_target_columns[[i]])$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column4"){
        print(xtabs(~experimental_column4,data=attributes(extended_target_columns[[i]])$newData))
        cat("\n")
      }
      if(target_columns_renamed[i]=="experimental_column5"){
        print(xtabs(~experimental_column5,data=attributes(extended_target_columns[[i]])$newData))
        cat("\n")
      }

    }

    ###### power curve simulation
    for(i in 1:length(target_columns)){

      if(levels[i]==1){
        pc=simr::powerCurve(extended_target_columns[[i]], test=simr::fixed("condition_column"),along=target_columns_renamed[i], nsim=nsimn,
                            breaks=seq(1,max_size[i],breaks[i])   )
      }else{
        pc=simr::powerCurve(extended_target_columns[[i]], test=simr::fixed("condition_column"),within=target_columns_renamed[i], nsim=nsimn,
                            breaks=seq(1,max_size[i],breaks[i])   )
      }





      if(length(output)==0){
        jpeg(paste("Power_curve_experiment_",target_columns[i],"_",Sys.Date(),".jpeg",sep=""))
      }else{
        jpeg(output[i])
      }

      print(plot(pc))
      print(mtext(response_column))
      dev.off()
      plot(pc)

    }

  }


}

