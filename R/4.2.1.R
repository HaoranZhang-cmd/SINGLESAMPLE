#' Draw an empirical ROC curve and calculate AUC for ordinal data
#' Also compute variance and CI for AUC
#' @param data a matrix with two rows.
#'             The first row is test results for diseased subjects
#'             The second row is test results for nondiseased subjects
#' @return an empirical ROC curve
#' #example
#' #data<-matrix(c(1,2,3,14,42,38,25,15,19,4),nrow=2,byrow=TRUE)
#' #roc.ordinal.emp(data)
#' @export
#-----4.2.1 Draw a scatter plot
roc.ordinal.emp<-function(data){
  if(dim(data)[1]!=2){
    stop("The row for table data should be 2")
  }
  K<-dim(data)[2]
  n1<-sum(data[1,])
  n0<-sum(data[2,])
  se<-c(rep(NA, length = K),0)
  FPR<-c(rep(NA, length = K),0)
  for(i in 1:K){
    se[i]= sum(data[1,i:K])/n1
    FPR[i] = sum(data[2,i:K])/n0
  }
  plot(FPR,se, type="b", xlab="FPR", ylab="Sensitivity", col="green",main=paste("Empirical ROC curve for Ordinal Data with ", K, "Categories"))
}
