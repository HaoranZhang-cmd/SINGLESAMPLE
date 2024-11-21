#' Function to calculate the full area and partial area and their variances
#' @param data0 ordinal test results of nondiseased subjects
#' @param data1 ordinal test results of diseased subjects
#' @param e1 lower bound of FPR for partial area
#' @param e2 upper bound of FPR for partial area
#' @return estimation of full area,partial area and their variances
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #roc.CI.ordinal(data0,data1,0,0.2)
#' @export
roc.CI.ordinal<-function(data0,data1,e1,e2){
  estimation<-params_estimate(data0,data1)
  a<-estimation$a
  b<-estimation$b
  var_a<-estimation$cov_matrix[1,1]
  var_b<-estimation$cov_matrix[2,2]
  covar_ab<-estimation$cov_matrix[1,2]
  #-------------Area and Partial Area under the ROC Curve (Parametric Methods)
  area_func<-function(x){
    return(pnorm(a+b*qnorm(x)))
  }
  area.full<-integrate(area_func,0,1)$value
  area.partial<-integrate(area_func,e1,e2)$value
  h1<-(qnorm(e1)+a*b/(1+b^2))*sqrt(1+b^2)
  h2<-(qnorm(e2)+a*b/(1+b^2))*sqrt(1+b^2)
  f<-exp(-a^2/2/(1+b^2))/sqrt(2*pi*(1+b^2))
  g<--a*b*exp(-a^2/2/(1+b^2))/sqrt(2*pi*(1+b^2)^3)
  partialf<-exp(-a^2/2/(1+b^2))/sqrt(2*pi*(1+b^2))*(pnorm(h2)-pnorm(h1))
  partialg<-exp(-a^2/2/(1+b^2))*(exp(-h1^2/2)-exp(-h2^2/2))/(2*pi*(1+b^2))-a*b*exp(-a^2/2/(1+b^2))*(pnorm(h2)-pnorm(h1))/(sqrt(2*pi*(1+b^2)^3))
  var_full<-f^2*var_a+g^2*var_b+2*f*g*covar_ab
  var_partial<-partialf^2*var_a+partialg^2*var_b+2*partialf*partialg*covar_ab
  output<-list(area.full=NULL,var_full=NULL,area.partial=NULL,var_partial=NULL)
  output$area.full<-area.full
  output$var_full<-var_full
  output$area.partial<-area.partial
  output$var_partial<-var_partial
  return(output)
}

#' function to estimate confidence intervals for full area and partial area
#' @param data0 ordinal test results of nondiseased subjects
#' @param data1 ordinal test results of diseased subjects
#' @param e1 lower bound of FPR for partial area
#' @param e2 upper bound of FPR for partial area
#' @param alpha significance level
#' @return CI_full,CI_partial:directly using Wald's interval
#'         CI_full_logit,CI_partial_logit:using logit transformation together with Wald's interval
#'         CI_partial_McClish:using McClish's transformation
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #CI_estimation_ordinal(data0,data1,0,0.2,0.05)
#' @export
CI_estimation_ordinal<-function(data0,data1,e1,e2,alpha){
  logit<-function(x){
    return(log(x/(1-x)))
  }
  result<-roc.CI.ordinal(data0,data1,e1,e2)
  A<-result$area.full
  varA<-result$var_full
  CI_full<-c(A-qnorm(1-alpha/2)*sqrt(varA),A+qnorm(1-alpha/2)*sqrt(varA))
  LL<-logit(A)-qnorm(1-alpha/2)*sqrt(varA)/(A*(1-A))
  UL<-logit(A)+qnorm(1-alpha/2)*sqrt(varA)/(A*(1-A))
  CI_full_logit<-c(exp(LL)/(1+exp(LL)),exp(UL)/(1+exp(UL)))

  #using the partial area index
  Amx<-e2-e1
  A<-result$area.partial
  varA<-result$var_partial
  CI_partial<-c(A-qnorm(1-alpha/2)*sqrt(varA),A+qnorm(1-alpha/2)*sqrt(varA))
  varlogit<-Amx^2*varA/(A*(Amx-A))^2
  LL<-logit(A/Amx)-qnorm(1-alpha/2)*sqrt(varlogit)
  UL<-logit(A/Amx)+qnorm(1-alpha/2)*sqrt(varlogit)
  CI_partial_logit<-c(Amx*exp(LL)/(1+exp(LL)),Amx*exp(UL)/(1+exp(UL)))
  #McClish transformation
  varMcClish<-4*Amx^2*varA/(Amx^2-A^2)^2
  LL<-log((Amx+A)/(Amx-A))-qnorm(1-alpha/2)*sqrt(varMcClish)
  UL<-log((Amx+A)/(Amx-A))+qnorm(1-alpha/2)*sqrt(varMcClish)
  CI_partial_McClish<-c(Amx*(exp(LL)-1)/(exp(LL)+1),Amx*(exp(UL)-1)/(exp(UL)+1))
  output<-list(CI_full=NULL,CI_partial=NULL,CI_full_logit=NULL,CI_partial_logit=NULL,CI_partial_McClish=NULL)
  output$CI_full<-CI_full
  output$CI_partial<-CI_partial
  output$CI_full_logit<-CI_full_logit
  output$CI_partial_logit<-CI_partial_logit
  output$CI_partial_McClish<-CI_partial_McClish
  return(output)
}

