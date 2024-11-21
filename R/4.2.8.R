#' function to estimate the parameters correlated to bigamma model
#' @param data0 ordinal test results of nondiseased subjects
#' @param data1 ordinal test results of diseased subjects
#' @return estimation of all thresholds,the bigamma parameters alpha and sigma for nondiseased subjects,the bigamma parameters alpha and sigma for diseased subjects
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #params_estimate_bigamma(data0,data1)
#' @export
params_estimate_bigamma<-function(data0,data1){
  loglikelihood<-function(params)
  {
    n0<-sum(data0)
    n1<-sum(data1)
    K<-length(data0)
    th<-params[1:(K-1)]
    alpha<-params[K]
    sigma0<-params[K+1]
    sigma1<-params[K+2]
    if(any(diff(th)<=0)) {
      return(Inf)
    }
    if(alpha<=0){
      return(Inf)
    }
    if(sigma0<=0||sigma1<=0){
      return(Inf)
    }
    result<-data0[1]*log(pgamma(th[1],shape=alpha,scale=sigma0),base=exp(1))+data0[K]*log(1-pgamma(th[K-1],shape=alpha,scale=sigma0),base=exp(1))+data1[1]*log(pgamma(th[1],shape=alpha,scale=sigma1),base=exp(1))+data1[K]*log(1-pgamma(th[K-1],shape=alpha,scale=sigma1),base=exp(1))
    for(i in 1:(K-2))
    {
      result<-result+data0[i+1]*log(pgamma(th[i+1],shape=alpha,scale=sigma0)-pgamma(th[i],shape=alpha,scale=sigma0), base=exp(1))+data1[i+1]*log(pgamma(th[i+1],shape=alpha,scale=sigma1)-pgamma(th[i],shape=alpha,scale=sigma1), base=exp(1))
    }
    return(-result)
  }
  result<-nlminb(start=c(seq(0.1,1.1,length=length(data0)-1),1,1,1),loglikelihood)
  print(result)
  estimation<-result$par
  return(estimation)
}

#' function to plot the ROC curve under bigamma assumption
#' @param data0 test results of nondiseased subjects
#' @param data1 test results of diseased subjects
#' @return a smooth ROC curve
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #roc.bigamma(data0,data1)
#' @export
roc.bigamma<-function(data0,data1){
  estimation<-params_estimate_bigamma(data0,data1)
  print(estimation)
  K<-length(data0)
  alpha<-estimation[K]
  sigma0<-estimation[K+1]
  sigma1<-estimation[K+2]
  c<-seq(0,1000,by=0.1)
  x<-rep(0,length=length(c))
  y<-rep(0,length=length(c))
  for(i in 1:length(c)){
    x[i]<-pgamma(c[i],shape=alpha,scale=sigma0,lower.tail=FALSE)
    y[i]<-pgamma(c[i],shape=alpha,scale=sigma1,lower.tail=FALSE)
  }
  plot(x,y,main=c("Estimation of Smooth ROC Curve with",length(data0),"Categories"),xlab="FPR",ylab="TPR",col="green",type="l")
}


#' function to estimate the parameters correlated to logistic model
#' @param data0 ordinal test results of nondiseased subjects
#' @param data1 ordinal test results of diseased subjects
#' @return estimation of all thresholds,the logistic parameters mu and gamma for nondiseased subjects,the logistic parameters mu and gamma for diseased subjects
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #params_estimate_logistic(data0,data1)
#' @export
params_estimate_logistic<-function(data0,data1){
  loglikelihood<-function(params)
  {
    n0<-sum(data0)
    n1<-sum(data1)
    K<-length(data0)
    th<-params[1:(K-1)]
    mu0<-0
    gamma0<-1
    mu1<-params[K]
    gamma1<-params[K+1]
    if(any(diff(th)<=0)) {
      return(Inf)
    }
    if(gamma0<=0||gamma1<=0){
      return(Inf)
    }
    result<-data0[1]*log(plogis(th[1],location=mu0,scale=gamma0),base=exp(1))+data0[K]*log(1-plogis(th[K-1],location=mu0,scale=gamma0),base=exp(1))+data1[1]*log(plogis(th[1],location=mu1,scale=gamma1),base=exp(1))+data1[K]*log(1-plogis(th[K-1],location=mu1,scale=gamma1),base=exp(1))
    for(i in 1:(K-2))
    {
      result<-result+data0[i+1]*log(plogis(th[i+1],location=mu0,scale=gamma0)-plogis(th[i],location=mu0,scale=gamma0), base=exp(1))+data1[i+1]*log(plogis(th[i+1],location=mu1,scale=gamma1)-plogis(th[i],location=mu1,scale=gamma1), base=exp(1))
    }
    return(-result)
  }
  result<-nlminb(start=c(seq(0,1,length=length(data0)-1),0,1),loglikelihood)
  print(result)
  estimation<-result$par
  return(estimation)
}

#' function to plot the ROC curve under logistic assumption
#' @param data0 test results of nondiseased subjects
#' @param data1 test results of diseased subjects
#' @return a smooth ROC curve
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #roc.logistic(data0,data1)
#' @export
roc.logistic<-function(data0,data1){
  estimation<-params_estimate_logistic(data0,data1)
  print(estimation)
  K<-length(data0)
  mu0<-0
  gamma0<-1
  mu1<-estimation[K]
  gamma1<-estimation[K+1]
  c<-seq(-1000,1000,by=0.1)
  x<-rep(0,length=length(c))
  y<-rep(0,length=length(c))
  for(i in 1:length(c)){
    x[i]<-plogis(c[i],location=mu0,scale=gamma0,lower.tail=FALSE)
    y[i]<-plogis(c[i],location=mu1,scale=gamma1,lower.tail=FALSE)
  }
  plot(x,y,main=c("Estimation of Smooth ROC Curve with",length(data0),"Categories"),xlab="FPR",ylab="TPR",col="green",type="l")
}
