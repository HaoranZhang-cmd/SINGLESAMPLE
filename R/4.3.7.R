#' FPR of the optimal decision point using the gaussian kernel
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param m it equals R*(1-p)/p,where p is prevalence rate and R is ratio of the cost
#' @return the optimal FPR and the corresponding TPR
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.gaussian_optimal(data0,data1,1)
#' @export
roc.gaussian_optimal<-function(data0,data1,m)
{
  data0<-roc.transformed(data0,data1)$data0_transformed
  data1<-roc.transformed(data0,data1)$data1_transformed
  x<-rep(0,length=2001)
  y<-rep(0,length=2001)
  z<-rep(0,length=2001)
  n0<-length(data0)
  h0<-0.9*min(sd(data0),(sort(data0)[ceiling(length(data0)*0.75)]-sort(data0)[ceiling(length(data0)*0.25)])/1.34)/(n0)^(1/5)
  for(i in -1000:1000)
  {
    result0<-0
    for(j in 1:n0)
    {
      result0<-result0+pnorm((data0[j]-i)/h0)/n0
    }
    x[i+1001]<-result0
  }
  n1<-length(data1)
  h1<-0.9*min(sd(data1),(sort(data1)[ceiling(length(data1)*0.75)]-sort(data1)[ceiling(length(data1)*0.25)])/1.34)/(n1)^(1/5)
  for(i in -1000:1000)
  {
    result1<-0
    for(j in 1:n1)
    {
      result1<-result1+pnorm((data1[j]-i)/h1)/n1
    }
    y[i+1001]<-result1
  }
  for(i in 1:2001){
    z[i]<-y[i]-m*x[i]
  }
  return(c(x[which.max(z)],y[which.max(z)]))
}

#' FPR of the optimal decision point using the biweight kernel
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param m it equals R*(1-p)/p,where p is prevalence rate and R is ratio of the cost
#' @return the optimal FPR and the corresponding TPR
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.biweight_optimal(data0,data1,1)
#' @export
roc.biweight_optimal<-function(data0,data1,m)
{
  data0<-roc.transformed(data0,data1)$data0_transformed
  data1<-roc.transformed(data0,data1)$data1_transformed
  n0<-length(data0)
  h0<-0.9*min(sd(data0),(sort(data0)[ceiling(length(data0)*0.75)]-sort(data0)[ceiling(length(data0)*0.25)])/1.34)/(n0)^(1/5)
  x<-rep(0,length=2001)
  y<-rep(0,length=2001)
  z<-rep(0,length=2001)
  for(i in -1000:1000)
  {
    result0<-0
    for(j in 1:n0)
    {
      if(i<data0[j]+h0)
      {
        f<-function(t)
        {
          return(15*(1-((t-data0[j])/h0)^2)^2/(16*n0*h0))
        }
        k<-integrate(f,lower=max(i,data0[j]-h0),upper=data0[j]+h0)
        result0<-result0+k$value
      }
    }
    x[i+1001]<-result0
  }

  n1<-length(data1)
  h1<-0.9*min(sd(data1),(sort(data1)[ceiling(length(data1)*0.75)]-sort(data1)[ceiling(length(data1)*0.25)])/1.34)/(n1)^(1/5)
  for(i in -1000:1000)
  {
    result1<-0
    for(j in 1:n1)
    {
      if(i<data1[j]+h1)
      {
        f<-function(t)
        {
          return(15*(1-((t-data1[j])/h1)^2)^2/(16*n1*h1))
        }
        k<-integrate(f,lower=max(i,data1[j]-h1),upper=data1[j]+h1)
        result1<-result1+k$value
      }
    }
    y[i+1001]<-result1
  }
  for(i in 1:2001){
    z[i]<-y[i]-m*x[i]
  }
  return(c(x[which.max(z)],y[which.max(z)]))
}

#' FPR of the optimal decision point when assuming binormality
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param m it equals R*(1-p)/p,where p is prevalence rate and R is ratio of the cost
#' @return the optimal FPR and the corresponding TPR and decision threshold
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.binormality_optimal(data0,data1,1)
#' @export
roc.binormality_optimal<-function(data0,data1,m){
  result<-list(binormality=NULL,threshold=NULL)
  a<-roc.continuous.smooth(data0,data1)$a
  b<-roc.continuous.smooth(data0,data1)$b
  if(b==1){
    FPR_optimal<-pnorm(-a/2-log(m,base=exp(1))/a)
    TPR_optimal<-pnorm(1/2-log(m,base=exp(1))/a)
  }
  else{
    FPR_optimal<-pnorm((a*b-sqrt(a^2+2*(1-b^2)*log(m/b)))/(1-b^2))
    TPR_optimal<-pnorm((a-b*sqrt(a^2+2*(1-b^2)*log(m/b)))/(1-b^2))
  }
  result$binormality<-c(FPR_optimal,TPR_optimal)
  result$threshold<-mean(data0)-sd(data0)*qnorm(FPR_optimal)
  return(result)
}
