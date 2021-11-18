#' @title rename_experimental_variables
#'
#' @description This function concatenate experimental variable names to reflect the nested experimental structure. For example, for Cell-lineA exists in ExperimentI and ExperimentII, the concatenated name will be ExperimentI_Cell-lineA.
#'
#'
#' @param data Input data
#' @param experimental_columns Name of the variable related to experimental design such as "experiment", "plate", and "cell_line".
#'
#'
#' @return Matrix with concatenated experimental variable names
#'
#' @export
#'
#' @examples rename_experimental_variables(data, c("experiment","line"))




rename_experimental_variables<-function(data, experimental_columns){



  print("Experimental column names should be sorted based on a hierarchy. ex) c('Experiment', 'plate', 'cell_line') ")

  i=length(experimental_columns)

  result=NULL
  if(i==2){

    result=as.matrix( paste(data[,experimental_columns[1]], data[,experimental_columns[2]], sep="_") )

  }else{
    while(i>1){
      result=cbind(.paste_columns(data[,experimental_columns[1:(i-1)]], as.matrix(data[,experimental_columns[i]]) ),result)
      i=i-1
    }

  }

  colnames(result)=paste0(experimental_columns[2:length(experimental_columns)],"_v2")
  return(cbind(result,data))


}



.paste_columns<-function(data1,data2){
  data1=as.matrix(data1)
  data2=as.matrix(data2)


  if( ncol(data1) == 1 & ncol(data2) == 1){
    if( sum(data1 == data2)==nrow(data1)){
      return(as.matrix(data1) )
    }
  }
  {
    return( as.matrix(paste(.paste_columns(data1[,1:(ncol(data1)-1)],data1[,ncol(data1)]),as.matrix(data2),sep="_") ) )
  }
}


