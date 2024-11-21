#' Compute Odds ratio and its confidence limits for binary-scale tests
#'
#' @param data 2*2 table of test results of diseased and nondiseased subjects
#' @param alpha significance level
#' @return estimation of odds ratio and its confidence limits
#' #example
#' #data<-matrix(c(22,3,2,3),nrow=2,byrow=TRUE)
#' #Oddsratio(data,0.05)
#' @export
Oddsratio<-function(data,alpha){
  if( dim(data)[1]!=2 & dim(data)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }
  if(data[1,1]*data[1,2]*data[2,1]*data[2,2]!=0){
    s1<-data[1,1]
    s0<-data[1,2]
    r1<-data[2,1]
    r0<-data[2,2]
  }
  else{
    s1<-data[1,1]+0.5
    s0<-data[1,2]+0.5
    r1<-data[2,1]+0.5
    r0<-data[2,2]+0.5
  }
  n1<-s1+s0
  n0<-r1+r0
  m1<-s1+r1
  m0<-s0+r0
  N<-n0+n1
  OR<-(s1*r0)/(s0*r1)
  lnOR.var<-1/s1+1/s0+1/r1+1/r0
  z<-qnorm(1-alpha/2)
  lowerlimit<-OR*exp(-z*sqrt(lnOR.var))
  upperlimit<-OR*exp(z*sqrt(lnOR.var))
  result<-list(OddsRatio=NULL,interval=NULL)
  result$OddsRatio<-OR
  result$interval<-rbind(lowerlimit,upperlimit)
  colnames(result$interval)<-'Estimation of confidence limits for Odds Ratio'
  rownames(result$interval)<-c('lower confidence limit','upper confidence limit')
  print(result)
}
