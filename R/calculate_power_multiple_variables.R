#' This function uses simulation to perform power analysis. It is designed to explore the power of biological experiments and to suggest an optimal number of experimental variables
#' with reasonable power. The backbone of the function is based on simr package, which fits a fixed effect or mixed effect model based on the observed data and simulates response
#' variables. Users can test the power of different combinations of experimental variables and parameters.
#'
#' Note: The current version does not accept categorical response variables, sample size parameters smaller than the observed samples size
#'
#' @import multtest
#' @import simr
#' @import lme4
#' @import lmerTest
#' @import readxl
#'
#' @param power_curve: 1: Power simulation over a range of sample sizes or levels. 0: Power calculation over a single sample size or a level (default).
#' @param data: Input data
#' @param variance_estimate_from: "data" (defalut) or "ICC". "data": If variances should be calculated from the input data. "ICC": User should provide Intra-Class Coefficients (ICC) in addition to observed data
#' @param condition_variable : Name of the condition varaible (ex variable with values such as control/case). The input file has to have a corresponding column name
#' @param experimental_variable : Name of the varaibles related to experimental desgin such as "experiment", "plate", and "cell_line".
#' @param response_variable : Name of the variable observed by performing the experiment. ex) intensity.
#' @param categorical_condition_variable: Specify whether the condition variable is categorical. 1: Categorical (default), 0: Continuous.
#' @param categorical_response_variable: Default: the observed variable is continous. Categorical response variable will be implemented in the future. 1: Categorical , 0: Continuous (default).
#' @param nsim: The number of simulations to run. Default=1000
#' @param  output: Output file name
##### If power_curve = 1
#' @param  target_parameter: Name of the experimental parameter to use for the power calculation.
#' @param  levels: 1: Amplify the number of corresponding target parameter. 0: Amplify the number of samples from the corresponding target parameter, ex) If target_parameter = c("experiment","cell_line") and if you want to expand the number of experiment and sample more cells from each cell line, levels = c(1,0).
#' @param  max_size: Maximum levels or sample sizes to test. Default: the current level or the current sample size x 5. ex) If max_levels = c(10,5), it will test upto 10 experiments and 5 cell lines.
#' @param  breaks: Levels /sample sizes of the variable to be specified along the power curve.. Default: max(1, round( the number of current levels / 5 ))
##### If variance_estimate_from = "data"
#' @param  effect_size: If you know the effect size of your condition variable, provided it. If the effect size is not provided, it will be estimated from your data
##### If variance_estimate_from = "varaince"
#' @param  ICC: Intra-Class Coefficients (ICC) for each parameter
#' @return A power curve image or a power calculation result printed in a text file
#'
#' @export
calculate_power_pdf <- function(power_curve=0,data,variance_estimate_from="data",condition_variable, experimental_variable, response_variable,
                            categorical_condition_variable=1, categorical_response_variable=0, nsimn=1000,
                            target_parameters, levels=NULL, max_size=NULL, breaks=NULL, effect_size=NULL, ICC=NULL, output=NULL){

  if(length(levels)!=length(target_parameters)){ print("User should specify levels of all target parameters") }
  #  variance_estimate_from="data"; levels=NULL; max_size=NULL; breaks=NULL; categorical_condition_variable=1; categorical_response_variable=1;nsimn=1000; effect_size=NULL; ICC=NULL
  #  data=inputData;variance_estimate_from="data";condition_variable="classification";experimental_variable=c("experiment","plate","line"); response_variable="perimeter";target_parameters="experiment"; levels=1; categorical_response_variable=0; nsimn=100
  # #
  #  power_curve=1;data=inputData2
  #  variance_estimate_from="ICC";condition_variable="classification";
  #  experimental_variable=c("experiment","exp_plate","line"); response_variable="perimeter"
  #  nsimn=10; target_parameters=c("experiment"); levels=1; ICC=c(0.2,0.15,0.3)
  # # # load.libs <- c(
  #   "multtest",
  #   "simr",
  #   "lmerTest",
  #   "readxl")
  # options(install.packages.check.source = "no")
  # options(install.packages.compile.from.source = "never")
  # if (!require("pacman")) install.packages("pacman"); library(pacman)
  # p_load(load.libs, update = TRUE, character.only = TRUE)
  # status <- sapply(load.libs,require,character.only = TRUE)
  # if(all(status)){
  #   print("SUCCESS: You have successfully installed and loaded all required libraries.")
  # } else{
  #   cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  #   status
  # }
  #
  require(multtest)
  require(simr)
  require(lme4)
  require(lmerTest)
  require(readxl)


  ####### read data

  # if(grepl("\\.csv",data)){
  #   Data=read.csv(data)
  # }else if(grepl("\\.xlsx",data)){
  #   Data=readxl::read_xlsx(data)
  # }else if(grepl("\\.xls$",data)){
  #   Data=readxl::read_xls(data)
  # }else if(grepl("\\.txt",data)){
  #   Data=read.table(data,header=TRUE)
  # }



  ####### remove empty lines
  Data <- data[complete.cases(data),]

  cat("\n")
  print("__________________________________________________________________Summary of data:")
  print(summary(Data))
  cat("\n")

  colnames_original=colnames(Data)
  experimental_variable_index=NULL
  ####### assign categorical variables
  if(categorical_condition_variable==1) Data[,condition_variable]=as.factor(Data[,condition_variable])
  if(categorical_response_variable==1) {
    cat("\n")
    print("_________________________________Categorical response variable is not accepted in the current version")}#Data[,response_variable]=as.factor(Data[,response_variable])
  cat("\n")
  for(i in 1:length(experimental_variable)){
    Data[,experimental_variable[i]]=as.factor(Data[,experimental_variable[i]])
    experimental_variable_index=c(experimental_variable_index,which(colnames(Data)==experimental_variable[i]))
    colnames(Data)[experimental_variable_index[i]]=paste("experimental_variable",i,sep="")
    cat("\n")
    print(paste("_________________________________",experimental_variable[i]," is assigned to experimental_variable",i,sep=""))
    cat("\n")
  }
  colnames(Data)[which(colnames(Data)==condition_variable)]="condition_variable"
  colnames(Data)[which(colnames(Data)==response_variable)]="response_variable"

  ###### indices of target parameters in experimental variables
  target_i=NULL
  target_parameters_renamed=NULL
  ###### match target parameters
  for(i in 1:length(target_parameters)){
    cn=colnames(Data)[which(colnames_original==target_parameters[i])]
    target_parameters_renamed[i]=cn

    target_i=c(target_i,as.integer(substr(cn,nchar(cn),nchar(cn))))
  }

  ####### run the formula
  if(variance_estimate_from=="data"){
    if(length(experimental_variable)==1){
      lmerFit <- lme4::lmer(response_variable ~ condition_variable + (1 | experimental_variable1), data=Data)
    }else if(length(experimental_variable)==2){
      lmerFit <- lme4::lmer(response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2), data=Data)
    }else if(length(experimental_variable)==3){
      lmerFit <- lme4::lmer(response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2) + (1 | experimental_variable3), data=Data)
    }else if(length(experimental_variable)==4){
      lmerFit <- lme4::lmer(response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2) + (1 | experimental_variable3) + (1 | experimental_variable4), data=Data)
    }else if(length(experimental_variable)==5){
      lmerFit <- lme4::lmer(response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2) + (1 | experimental_variable3) + (1 | experimental_variable4) + (1 | experimental_variable5), data=Data)
    }
  }else{


    lmerFit=stats::lm(response_variable ~ condition_variable,data=Data)
  }





  slmerFit <- summary(lmerFit)
  cat("\n")
  print("__________________________________________________________________Model statistics:")
  print(slmerFit)
  cat("\n")

  #fixed_effects=c(lme4::fixef(lmerFit)["(Intercept)"],lme4::fixef(lmerFit)["condition_variable1"])
  fixed_effects=slmerFit$coefficients[,1]

  ##### ICC based variance estimation
  if(variance_estimate_from=="ICC"){

    ####### If there is only one category
    for(i in 1:length(experimental_variable)){
      if(length(table(Data[,experimental_variable_index[i]]))==1){
        categroy_temp=paste0(Data[,experimental_variable_index[i]][1],"_2")
        levels(Data[,experimental_variable_index[i]])=c(levels(Data[,experimental_variable_index[i]]),categroy_temp)
        Data[,experimental_variable_index[i]][sample(1:length(Data[,1]),round(length(Data[,1])/2))]=rep(categroy_temp,round(length(Data[,1])/2))
      }
    }





    if(length(experimental_variable)==1){
      varEs=(slmerFit$sigma)^2*ICC/(1-ICC)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variable:",varEs))
      cat("\n")
    }else if(length(experimental_variable)==2){

      a <- matrix(c( 1-ICC[1], -ICC[2],-ICC[1],1-ICC[2]), nrow=2, ncol=2)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2 ) , nrow=2, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      artificial_glmer=simr::makeLmer(formula = response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2)
                                      , data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )

    }else if(length(experimental_variable)==3){

      a <- matrix(c( 1-ICC[1], -ICC[2], -ICC[3], -ICC[1], 1-ICC[2], -ICC[3], -ICC[1], -ICC[2], 1-ICC[3]), nrow=3, ncol=3)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2, ICC[3]*(slmerFit$sigma)^2 ) , nrow=3, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      #vars=list(experimental_variable1=varEs[1],experimental_variable2=varEs[2],experimental_variable3=varEs[3])
      #Data_temp=Data
      #Data=Data[,-which(colnames(Data)=="experimental_variable1")]
      #Data=Data[,-which(colnames(Data)=="experimental_variable2")]
      #Data=Data[,-which(colnames(Data)=="experimental_variable3")]
      artificial_glmer=simr::makeLmer(formula = response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2) +
                                        (1 | experimental_variable3) , data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )


    }else if(length(experimental_variable)==4){

      a <- matrix(c( 1-ICC[1], -ICC[2], -ICC[3], -ICC[4], -ICC[1], 1-ICC[2], -ICC[3], -ICC[4], -ICC[1], -ICC[2], 1-ICC[3], -ICC[4], -ICC[1], -ICC[2], -ICC[3], 1-ICC[4]), nrow=4, ncol=4)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2, ICC[3]*(slmerFit$sigma)^2, ICC[4]*(slmerFit$sigma)^2 ) , nrow=4, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      artificial_glmer=simr::makeLmer(formula = response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2) +
                                        (1 | experimental_variable3) + (1 | experimental_variable4) , data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )

    }else if(length(experimental_variable)==5){

      a <- matrix(c( 1-ICC[1], -ICC[2], -ICC[3], -ICC[4], -ICC[5], -ICC[1], 1-ICC[2], -ICC[3], -ICC[4], -ICC[5], -ICC[1], -ICC[2], 1-ICC[3], -ICC[4], -ICC[5]
                     , -ICC[1], -ICC[2], -ICC[3], 1-ICC[4], -ICC[5], -ICC[1], -ICC[2], -ICC[3], -ICC[4], 1-ICC[5]), nrow=4, ncol=4)
      b <- matrix(c( ICC[1]*(slmerFit$sigma)^2, ICC[2]*(slmerFit$sigma)^2, ICC[3]*(slmerFit$sigma)^2, ICC[4]*(slmerFit$sigma)^2, ICC[5]*(slmerFit$sigma)^2 ) , nrow=4, ncol=1)
      varEs=solve(a,b)
      cat("\n")
      print(paste("__________________________________________________________________Estimated variance of the experimental variables:",paste(varEs)))
      cat("\n")

      artificial_glmer=simr::makeLmer(formula = response_variable ~ condition_variable + (1 | experimental_variable1) + (1 | experimental_variable2) +
                                        (1 | experimental_variable3) + (1 | experimental_variable4) + (1 | experimental_variable5), data=Data,
                                      VarCorr = as.list(varEs), sigma = slmerFit$sigma,
                                      fixef=fixed_effects )

    }


    cat("\n")
    print(artificial_glmer)
    print("__________________________________________________________________Model statistics:")
    print(summary(artificial_glmer))
    cat("\n")
    lmerFit=artificial_glmer
  }else{##not ACC



  }




  ####### print observed levels and sample sizes
  maxs=NULL
  mins=NULL
  lens=NULL
  for(i in 1:length(experimental_variable)){
    cat("\n")
    print(paste("__________________________________________________________________Levels and sample sizes of",experimental_variable[i]))
    cat("\n")
    xtabs_s=stats::xtabs(~Data[,paste0("experimental_variable",i)])


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
    for(i in 1:length(condition_variable)){
      if(effect_size[i]!=-1)
        fixef(lmerFit)["condition_variable1"] <- effect_size[i]
      cat("\n")
      print(paste0("_________________________________Effect size of the condition_variable is now ",effect_size[i]))
      cat("\n")
    }
  }


  if(power_curve==0){

    ####### extend parameter levels and sample sizes

    if(length(max_size)==0){max_size=rep(0,length(target_parameters))}
    if(length(breaks)==0){
      breaks=rep(0,length(target_parameters))
    }

    for(i in 1:length(target_parameters)){
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
          extended_target_parameters=simr::extend(lmerFit,along=target_parameters_renamed[i],n=max_size[i])
        }else{
          extended_target_parameters=simr::extend(extended_target_parameters,along=target_parameters_renamed[i],n=max_size[i])
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
          extended_target_parameters=simr::extend(lmerFit,within=target_parameters_renamed[i],n=max_size[i])
        }else{
          extended_target_parameters=simr::extend(extended_target_parameters,within=target_parameters_renamed[i],n=max_size[i])
        }
      }
      cat("\n")
      print(paste("Input max size of",target_parameters[i]))
      print(max_size[i])
      cat("\n")

    }

    cat("\n")
    print("__________________________________________________________________Extended parameters:")
    for(i in 1:length(target_parameters)){
      if(target_parameters_renamed[i]=="experimental_variable1"){
        print(xtabs(~experimental_variable1,data=attributes(extended_target_parameters)$newData))
        #print(xtabs(~experimental_variable2,data=attributes(extended_target_parameters)$newData))
        #print(xtabs(~experimental_variable3,data=attributes(extended_target_parameters)$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable2"){
        #print(xtabs(~experimental_variable1,data=attributes(extended_target_parameters)$newData))
        print(xtabs(~experimental_variable2,data=attributes(extended_target_parameters)$newData))
        # print(xtabs(~experimental_variable3,data=attributes(extended_target_parameters)$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable3"){
        # print(xtabs(~experimental_variable1,data=attributes(extended_target_parameters)$newData))
        # print(xtabs(~experimental_variable2,data=attributes(extended_target_parameters)$newData))
        print(xtabs(~experimental_variable3,data=attributes(extended_target_parameters)$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable4"){
        print(xtabs(~experimental_variable4,data=attributes(extended_target_parameters)$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable5"){
        print(xtabs(~experimental_variable5,data=attributes(extended_target_parameters)$newData))
        cat("\n")
      }

    }
    ###### power simulation

    ps=simr::powerSim(extended_target_parameters, test=fixed("condition_variable"),nsim=nsimn)
    cat("\n")
    print("__________________________________________________________________Power simulation result:")
    print(ps)
    cat("\n")
    if(length(output)!=0){
      sink(output)
    }else{
      sink(paste("Power_curve_experiment_",target_parameters[i],"_",Sys.Date(),".txt",sep=""))
    }

    print(ps)
    sink()





  }else{#power curve =1


    ####### extend parameter levels and sample sizes
    extended_target_parameters=list()
    if(length(max_size)==0){max_size=rep(0,length(target_parameters))}
    if(length(breaks)==0){
      breaks=rep(0,length(target_parameters))
    }

    for(i in 1:length(target_parameters)){
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
        extended_target_parameters=c(extended_target_parameters,list(simr::extend(lmerFit,along=target_parameters_renamed[i],n=max_size[i])))
      }else{ # increase sample sizes
        if(max_size[i]==0){
          max_size[i]=maxs[target_i[i]]*5
        }
        if(breaks[i]==0) breaks[i] = max(1, round( maxs[target_i[i]] / 5 ))
        extended_target_parameters=c(extended_target_parameters,list(simr::extend(lmerFit,within=target_parameters_renamed[i],n=max_size[i])))
      }
      cat("\n")
      print(paste("_________________________________Power simulation will be performed based on the max size of",target_parameters[i],":"))
      print(max_size[i])
      cat("\n")

    }


    cat("\n")
    print("__________________________________________________________________Extended parameters:")
    for(i in 1:length(target_parameters)){
      if(target_parameters_renamed[i]=="experimental_variable1"){
        print(xtabs(~experimental_variable1,data=attributes(extended_target_parameters[[i]])$newData))
        #print(xtabs(~experimental_variable2,data=attributes(extended_target_parameters[[i]])$newData))
        #print(xtabs(~experimental_variable3,data=attributes(extended_target_parameters[[i]])$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable2"){
        #print(xtabs(~experimental_variable1,data=attributes(extended_target_parameters[[i]])$newData))
        print(xtabs(~experimental_variable2,data=attributes(extended_target_parameters[[i]])$newData))
        # print(xtabs(~experimental_variable3,data=attributes(extended_target_parameters[[i]])$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable3"){
        # print(xtabs(~experimental_variable1,data=attributes(extended_target_parameters[[i]])$newData))
        # print(xtabs(~experimental_variable2,data=attributes(extended_target_parameters[[i]])$newData))
        print(xtabs(~experimental_variable3,data=attributes(extended_target_parameters[[i]])$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable4"){
        print(xtabs(~experimental_variable4,data=attributes(extended_target_parameters[[i]])$newData))
        cat("\n")
      }
      if(target_parameters_renamed[i]=="experimental_variable5"){
        print(xtabs(~experimental_variable5,data=attributes(extended_target_parameters[[i]])$newData))
        cat("\n")
      }

    }

    ###### power curve simulation
    for(i in 1:length(target_parameters)){

      if(levels[i]==1){
        pc=simr::powerCurve(extended_target_parameters[[i]], test=fixed("condition_variable"),along=target_parameters_renamed[i], nsim=nsimn,
                            breaks=seq(1,max_size[i],breaks[i])   )
      }else{
        pc=simr::powerCurve(extended_target_parameters[[i]], test=fixed("condition_variable"),within=target_parameters_renamed[i], nsim=nsimn,
                            breaks=seq(1,max_size[i],breaks[i])   )
      }

      #print(pc)
    #  if(length(output)==0){
    #    jpeg(paste("Power_curve_experiment_",target_parameters[i],"_",Sys.Date(),".jpeg",sep=""))
    #  }else{
    #    jpeg(output[i])
    #  }
      print(plot(pc))
      print(mtext(response_variable))

    #  dev.off()
    #  plot(pc)
    }

  }




}


calculate_power_multiple_variables <- function(data,
                                                variance_estimate_from="data",condition_variable="classification",
                                                experimental_variable=c("experiment","plate","line"), response_variable,
                                                nsimn=10, target_parameters="experiment", levels=1, max_size=15,pdf_name="test.pdf"){


  pdf(pdf_name)

  for(i in 1:length(response_variable)){

    result=calculate_power_pdf(power_curve=1,data=data,
                           variance_estimate_from=variance_estimate_from,condition_variable=condition_variable,
                           experimental_variable=experimental_variable, response_variable=response_variable[i],
                           nsimn=nsimn, target_parameters=target_parameters, levels=levels, max_size=max_size)
  }

  dev.off()

}

