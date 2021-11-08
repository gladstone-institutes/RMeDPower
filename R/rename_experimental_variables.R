#' @title rename_experimental_variables
#'
#' @export



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

