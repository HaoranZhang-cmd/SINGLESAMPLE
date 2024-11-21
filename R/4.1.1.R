#' Compute Sensitivity and Specificity for binary scale tests
#' @param data 2*2 table of test results of diseased and nondiseased subjects
#' @param alpha significance level
#' @return estimation of sensitivity and specificity,their variances,and the estimated
#' confidence intervals,including wald's,AC and ZL intervals.
#' #example
#' #data<-matrix(c(22,3,2,3),nrow=2,byrow=TRUE)
#' #roc.binary(data,0.05)
#' @export
#4.1.1 Sensitivity and Specificity
roc.binary<-function(data, alpha){
  if( dim(data)[1]!=2 || dim(data)[2]!=2 ){
    stop("The row and column for table data should be 2")
  }

  s1<-data[1,1]
  s0<-data[1,2]
  r1<-data[2,1]
  r0<-data[2,2]
  K<-dim(data)[2]
  n1<-sum(data[1,])
  n0<-sum(data[2,])

  result<-list(Estimate=NULL, Variance=NULL, Interval=NULL)

  #-----estimate ---------
  se<-s1/n1
  sp<-r0/n0
  #-----estimate variance ----
  se.var<-s1*s0/n1^3
  sp.var<-r1*r0/n0^3

  #-----interval ----
  #1.THE ORDINAL METHOD FOR CI
  #-----wald's interval ---------
  z<-qnorm(1-alpha/2)
  se.CI_wald<-c(se-z*sqrt(se.var), se+z*sqrt(se.var))
  sp.CI_wald<-c(sp-z*sqrt(sp.var), sp+z*sqrt(sp.var))

  #2.WHEN SAMPLESIZE IS LOW OR AUCCURANCY CLOSE 0/1 FOR ADJUSTED CI
  #-----AC interval ------------
  se.lower<-(s1+z^2/2)/(n1+z^2)-z*sqrt(s1*s0/n1+z^2/4)/(n1+z^2)
  se.upper<-(s1+z^2/2)/(n1+z^2)+z*sqrt(s1*s0/n1+z^2/4)/(n1+z^2)
  sp.lower<-(r0+z^2/2)/(n0+z^2)-z*sqrt(r0*r1/n0+z^2/4)/(n0+z^2)
  sp.upper<-(r0+z^2/2)/(n0+z^2)+z*sqrt(r0*r1/n0+z^2/4)/(n0+z^2)
  se.CI_AC<-c(se.lower,se.upper)
  sp.CI_AC<-c(sp.lower,sp.upper)

  zl<-qnorm(1-alpha/2)
  zu<-qnorm(alpha/2)
  n<-n0+n1
  loglink_se<-function(x){
    gamma<-(1-2*se)/sqrt(se*(1-se))
    a<-(-1)*sqrt(n)/(gamma/6)
    b<-gamma/2*(x/sqrt(n)-gamma/(6*n))
    return(a*((1-b)^(1/3)-1))
  }
  loglink_sp<-function(x){
    gamma<-(1-2*sp)/sqrt(sp*(1-sp))
    a<-(-1)*sqrt(n)/(gamma/6)
    b<-gamma/2*(x/sqrt(n)-gamma/(6*n))
    return(a*((1-b)^(1/3)-1))
  }
  se.CI_ZL_trans<-c(log(se/(1-se))-n^(-1/2)*(se*(1-se))^(-1/2)*loglink_se(zl),log(se/(1-se))-n^(-1/2)*(se*(1-se))^(-1/2)*loglink_se(zu))
  sp.CI_ZL_trans<-c(log(sp/(1-sp))-n^(-1/2)*(sp*(1-sp))^(-1/2)*loglink_sp(zl),log(sp/(1-sp))-n^(-1/2)*(sp*(1-sp))^(-1/2)*loglink_sp(zu))
  se.CI_ZL<-exp(se.CI_ZL_trans)/(1+exp(se.CI_ZL_trans))
  sp.CI_ZL<-exp(sp.CI_ZL_trans)/(1+exp(sp.CI_ZL_trans))

  result$Estimate<-rbind(se, sp)
  result$Variance<-rbind(se.var, sp.var)
  result$Interval<-rbind(se.CI_wald,se.CI_AC,se.CI_ZL,sp.CI_wald,sp.CI_AC,sp.CI_ZL)
  colnames(result$Estimate)[1]<-'estimate'
  colnames(result$Variance)[1]<-'variance'
  colnames(result$Interval)<-c('low_level','up_level')

  return(result)
}
