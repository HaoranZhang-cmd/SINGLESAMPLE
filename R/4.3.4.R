#--------------4.3.4 Area and Partial Area Under the ROC Curve
#' this function is used to calculate AUC for continuous data.It is not important.
#' @param x not important
#' @param y not important
#' @export
Psi<-function(x,y){
  if(x>y){
    return(1)
  }
  if(x==y){
    return(1/2)
  }
  if(x<y){
    return(0)
  }
}

#' function used to calculate AUC and partial AUC for continuous data using nonparametric methods
#' @param data0 continuous test result for undiseased subjects
#' @param data1 continuous test result for diseased subjects
#' @param e1 lower limit of FPR for partial AUC
#' @param e2 upper limit of FPR for partial AUC
#' @return ustimation of AUC,partial AUC and its variance using Dodd's and He's methods
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.area.nonparametric(data0,data1,0,0.2)
#' @export
roc.area.nonparametric<-function(data0,data1,e1,e2){
  result<-list(full_area=NULL,partial_area_Dodd=NULL,partial_area_He=NULL,var_partial_area=NULL)
  n0<-length(data0)
  h0<-0.9*min(sd(data0),(quantile(data0,probs=0.75)-quantile(data0,probs=0.25))/1.34)/(n0)^(1/5)
  n1<-length(data1)
  h1<-0.9*min(sd(data1),(quantile(data1,probs=0.75)-quantile(data1,probs=0.25))/1.34)/(n1)^(1/5)
  #--------------full area estimation
  full<-0
  for(i in 1:n1){
    for(j in 1:n0){
      full<-full+pnorm((data1[i]-data0[j])/sqrt(h1^2+h0^2))/(n0*n1)
    }
  }
  #--------------partial area estimation
  partial_Dodd<-0
  for(i in 1:n1){
    for(j in 1:n0){
      if(data1[i]>=data0[j]){
        if(data0[j]>quantile(data0,probs=e1)&&data0[j]<=quantile(data0,probs=e2)){
          partial_Dodd<-partial_Dodd+1/(n0*n1)
        }
      }
    }
  }
  partial_He<-0
  for(i in 1:n1){
    for(j in 1:n0){
      if(data0[j]>=quantile(data0,probs=e1)&&data0[j]<=quantile(data0,probs=e2)){
        partial_He<-partial_He+Psi(data1[i],data0[j])/(n0*n1)
      }
    }
  }
  #----------------variance estimation
  np<-0
  for(j in 1:n0){
    if(data0[j]>=quantile(data0,probs=e1)&&data0[j]<=quantile(data0,probs=e2)){
      np<-np+1
    }
  }
  tau<-np*partial_He/n0
  n<-n0+n1
  tmp1<-0
  tmp2<-0
  for(i in 1:n0){
    if(data0[i]>=quantile(data0,probs=e1)&&data0[i]<=quantile(data0,probs=e2)){
      V10<-0
      for(j in 1:n1){
        V10<-V10+Psi(data1[j],data0[i])/n0
      }
      tmp1<-tmp1+(V10-tau)^2
    }
  }
  for(j in 1:n1){
    V01<-0
    for(i in 1:n0){
      if(data0[i]>=quantile(data0,probs=e1)&&data0[i]<=quantile(data0,probs=e2)){
        V01<-V01+Psi(data1[j],data0[i])/np
      }
    }
    tmp2<-tmp2+(V01-tau)^2
  }
  var.He<-np*tmp1/(n*n*(np-1))+tmp2/(n1*(n1-1))
  #-----------------print results
  result$full_area<-full
  result$partial_area_Dodd<-partial_Dodd
  result$partial_area_He<-partial_He
  result$var_partial_area<-var.He
  print(result)
}


#' Function to calculate the full area and partial area and their variances under binormal assmuption
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param e1 lower bound of FPR for partial area
#' @param e2 upper bound of FPR for partial area
#' @return estimation of full area and partial area and their variances
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.var.continuous(data0,data1,0,0.2)
#' @export
roc.var.continuous<-function(data0,data1,e1,e2){
  estimation<-roc.continuous.smooth(data0,data1)
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


#' function to estimate confidence intervals for full area and partial area under binormal assumption
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @param e1 lower bound of FPR for partial area
#' @param e2 upper bound of FPR for partial area
#' @param alpha significance level
#' @return CI_full,CI_partial:directly using Wald's interval
#'         CI_full_logit,CI_partial_logit:using logit transformation together with Wald's interval
#'         CI_partial_McClish:using McClish's transformation
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #CI_estimation_continuous(data0,data1,0,0.2,0.05)
#' @export
CI_estimation_continuous<-function(data0,data1,e1,e2,alpha){
  logit<-function(x){
    return(log(x/(1-x)))
  }
  result<-roc.var.continuous(data0,data1,e1,e2)
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



