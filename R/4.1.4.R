#' Compute positive and negative LR as well as their confidence limits for binary-scale tests
#' @param data 2*2 table of test results of diseased and undiseased subjects
#' @param alpha significance level
#' @return Estimation of positive and negative LR,Estimation of the variance of logarithms of positive and negative LR,
#' Estimation of confidence limits for positive and negative LR
#' #example
#' #data<-matrix(c(22,3,2,3),nrow=2,byrow=TRUE)
#' #LRratio(data,0.05)
#' @export
LRratio<-function(data,alpha){
  if( dim(data)[1]!=2 & dim(data)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }
  s1<-data[1,1]
  s0<-data[1,2]
  r1<-data[2,1]
  r0<-data[2,2]
  n1<-sum(data[1,])
  n0<-sum(data[2,])
  m1<-sum(data[,1])
  m0<-sum(data[,2])
  N<-n0+n1

  Se<-s1/n1
  Sp<-r0/n0
  positiveLR<-Se/(1-Sp)
  lnpositiveLR.var<-(1-Se)/s1+Sp/r1
  z<-qnorm(1-alpha/2)
  positivelowerlimit<-positiveLR*exp(-z*sqrt(lnpositiveLR.var))
  positiveupperlimit<-positiveLR*exp(z*sqrt(lnpositiveLR.var))
  negativeLR<-(1-Se)/Sp
  lnnegativeLR.var<-Se/s0+(1-Sp)/r0
  negativelowerlimit<-negativeLR*exp(-z*sqrt(lnnegativeLR.var))
  negativeupperlimit<-negativeLR*exp(z*sqrt(lnnegativeLR.var))

  result<-list(Estimate1=NULL,Estimate2=NULL,Interval1=NULL,Interval2=NULL)
  result$Estimate1<-rbind(positiveLR,negativeLR)
  result$Estimate2<-rbind(lnpositiveLR.var,lnnegativeLR.var)
  result$Interval1<-rbind(positivelowerlimit,positiveupperlimit)
  result$Interval2<-rbind(negativelowerlimit,negativeupperlimit)

  colnames(result$Estimate1)[1]<-'Estimation of positive and negative LR'
  colnames(result$Estimate2)[1]<-'Estimation of the variance of logarithms of positive and negative LR'
  colnames(result$Interval1)[1]<-'Estimation of confidence limits for positive LR'
  colnames(result$Interval2)[1]<-'Estimation of confidence limits for negative LR'
  rownames(result$Estimate1)<-c('positive LR','negative LR')
  rownames(result$Estimate2)<-c('the variance of the logarithm of positive LR','the variance of the logarithm of positive LR')
  rownames(result$Interval1)<-c('lower confidence limit for positive LR','upper confidence limit for positive LR')
  rownames(result$Interval2)<-c('lower confidence limit for negative LR','upper confidence limit for negative LR')
  print(result)
}
