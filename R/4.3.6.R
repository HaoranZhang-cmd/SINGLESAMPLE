#-----------4.3.6 Fixed False Positive Rate - Sensitivity and the Decision Threshold
#' function to calculate sensitivity and decision threshold at a fixed FPR
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param SP the fixed specificity
#' @param alpha significance level
#' @return estimation of Se corresponding to the threshold at a fixed Sp,confidence intervals for estimated Se,
#' the decision threshold,the lower confidence limit for Se
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #SDT(data0,data1,0.9,0.05)
#' @export
SDT<-function(data0,data1,SP,alpha){
  n0<-length(data0)
  n1<-length(data1)
  TSP<-sort(data0)[ceiling(n0*SP)]
  Se.Greenhouse<-0
  for(i in 1:n1){
    if(data1[i]>TSP){
      Se.Greenhouse<-Se.Greenhouse+1/n1
    }
  }
  print(c("Estimation of sensitivity corresponding to this threshold is",Se.Greenhouse))
  varSe.Greenhouse<-Se.Greenhouse*(1-Se.Greenhouse)/n1
  new_data0<-roc.transformed(data0,data1)$data0_transformed
  new_data1<-roc.transformed(data0,data1)$data1_transformed
  new_TSP<-sort(new_data0)[ceiling(n0*SP)]
  h0<-0.9*min(sd(new_data0),(sort(new_data0)[ceiling(n0*0.75)]-sort(new_data0)[ceiling(n0*0.25)])/1.34)/(n0)^(1/5)
  h1<-0.9*min(sd(new_data1),(sort(new_data1)[ceiling(n1*0.75)]-sort(new_data1)[ceiling(n1*0.25)])/1.34)/(n1)^(1/5)
  f0<-0
  f1<-0
  for(i in 1:n0){
    f0<-f0+dnorm((new_TSP-new_data0[i])/h0)/(n0*h0)
  }
  for(i in 1:n1){
    f1<-f1+dnorm((new_TSP-new_data1[i])/h1)/(n1*h1)
  }
  varTSP<-SP*(1-SP)/(f0^2*n0)
  varSe.Linnet<-varSe.Greenhouse+varTSP*f1^2

  #------------------Zhou and Qin confidence interval
  z<-qnorm(1-alpha/2)
  B<-1000
  RSP<-rep(0,length=B)
  VSP<-0
  for(b in 1:B){
    x0<-sample(data0,n0,replace=TRUE)
    x1<-sample(data1,n1,replace=TRUE)
    for(i in 1:n1){
      if(x1[i]>=sort(data0)[ceiling(n0*SP)]){
        RSP[b]<-RSP[b]+1
      }
    }
    RSP[b]<-(RSP[b]+z^2/2)/(n1+z^2)
  }
  interval<-c(mean(RSP)-z*sqrt(var(RSP)),mean(RSP)+z*sqrt(var(RSP)))

  c_nonparametric<-new_TSP+qnorm(sqrt(1-alpha))*sqrt(varTSP)
  lambda<-roc.transformed(data0,data1)$lambda
  #用Box-Cox变换的逆变换换回原始数据
  threshold<-(lambda*c_nonparametric+1)^(1/lambda)
  Se_nonparametric<-0
  for(i in 1:n1){
    if(new_data1[i]>c_nonparametric){
      Se_nonparametric<-Se_nonparametric+1/n1
    }
  }
  lower_nonparametric<-Se_nonparametric-qnorm(sqrt(1-alpha))*sqrt(Se_nonparametric*(1-Se_nonparametric)/n1)
  print("Confidence intervals for the estimated sensitivity,Zhou and Qin")
  print(interval)
  print("The decision threshold")
  print(threshold)
  print("The lower confidence bound for sensitivity")
  print(lower_nonparametric)
}






