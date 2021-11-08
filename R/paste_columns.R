#' @title paste_columns
#'
#' @export


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
