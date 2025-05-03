#' Compute PPV and NPV for binary scale tests
#' prevalence rate=x/n
#' @param data 2*2 table of test results of diseased and nondiseased subjects
#' @param alpha significance level
#' @param n the number of subjects evaluated in a random sample from the relevant population
#' @param x the number of diseased subjects in n subjects
#' @return 1) Estimated PPV and NPV in the study population(where the original data comes from)
#' and their confidence intervals(based on Wald's interval)
#' 2) Estimated PPV and NPV in a relevant population(with a different prevalence rate as the study population)
#' and three different confidence intervals, the first is based on logit transformation,
#' the second is based on Delta method, and the third is based on objective Bayesian method.
#' #example
#' #data<-matrix(c(56,6,23,78),nrow=2,byrow=TRUE)
#' #roc.ppvnpv(data,0.05,124,632)
#' @export
#4.1.2 Predictive Value of a Positive or Negative
roc.ppvnpv<-function(data,alpha,x,n){
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
  #------------- Estimated value in the study population(where the original data comes from)
  Se<-s1/n1
  Sp<-r0/n0
  output<-list(Estimate_study_population=NULL,PPV.interval_study_population=NULL,
               NPV.interval_study_population=NULL,Estimate_relevant_population=NULL,
               PPV.interval1_relevant_population=NULL,PPV.interval2_relevant_population=NULL,
               PPV.interval3_relevant_population=NULL,NPV.interval1_relevant_population=NULL,
               NPV.interval2_relevant_population=NULL,NPV.interval3_relevant_population=NULL)
  p<-x/n
  output$Estimate_study_population<-c(s1/m1,r0/m0)
  var_study_population<-c(s1*r1/m1^3,r0*s0/m0^3)
  output$PPV.interval_study_population<-c(s1/m1-qnorm(1-alpha/2)*sqrt(s1*r1/m1^3),s1/m1+qnorm(1-alpha/2)*sqrt(s1*r1/m1^3))
  output$NPV.interval_study_population<-c(r0/m0-qnorm(1-alpha/2)*sqrt(r0*s0/m0^3),r0/m0+qnorm(1-alpha/2)*sqrt(r0*s0/m0^3))

  #Assume Se and Sp is the same, transform to a relevant population with different prevalence rate
  PPV<-p*Se/((p*Se)+(1-Sp)*(1-p))
  NPV<-Sp*(1-p)/(Sp*(1-p)+(1-Se)*p)
  #---------------three different confidence intervals
  #based on the logit transformation of the Bayes formula
  PPV.logit<-log(p*Se/((1-Sp)*(1-p)),base=exp(1))
  NPV.logit<-log((1-p)*Sp/((1-Se)*p),base=exp(1))
  PPV.varlogit<-(1-Se)/(Se*n1)+Sp/(n0*(1-Sp))
  NPV.varlogit<-Se/((1-Se)*n1)+(1-Sp)/(Sp*n0)
  z<-qnorm(1-alpha/2)
  PPV.interval1<-c(exp(PPV.logit-z*sqrt(PPV.varlogit))/(1+exp(PPV.logit-z*sqrt(PPV.varlogit))),exp(PPV.logit+z*sqrt(PPV.varlogit))/(1+exp(PPV.logit+z*sqrt(PPV.varlogit))))
  NPV.interval1<-c(exp(NPV.logit-z*sqrt(NPV.varlogit))/(1+exp(NPV.logit-z*sqrt(NPV.varlogit))),exp(NPV.logit+z*sqrt(NPV.varlogit))/(1+exp(NPV.logit+z*sqrt(NPV.varlogit))))
  #based on delta method
  p.var<-x*(n-x)/n^3
  Se.var<-s1*s0/n1^3
  Sp.var<-r1*r0/n0^3
  PPV.partialp<-Se*(1-Sp)/((p*Se+(1-p)*(1-Sp))^2)
  PPV.partialSe<-(1-p)*p*(1-Sp)/((p*Se+(1-p)*(1-Sp))^2)
  PPV.partialSp<-(p-1)*p*Se/((p*Se+(1-p)*(1-Sp))^2)
  PPV.var<-(PPV.partialp)^2*p.var+(PPV.partialSe)^2*Se.var+(PPV.partialSp)^2*Sp.var
  NPV.partialp<-(Se-1)*Sp/(Sp*(1-p)+(1-Se)*p)^2
  NPV.partialSe<-p*Sp*(1-p)/(Sp*(1-p)+(1-Se)*p)^2
  NPV.partialSp<-(1-Se)*p*(1-p)/(Sp*(1-p)+(1-Se)*p)^2
  NPV.var<-(NPV.partialp)^2*p.var+(NPV.partialSe)^2*Se.var+(NPV.partialSp)^2*Sp.var
  z<-qnorm(1-alpha/2)
  PPV.interval2<-c(PPV-z*sqrt(PPV.var),PPV+z*sqrt(PPV.var))
  NPV.interval2<-c(NPV-z*sqrt(NPV.var),NPV+z*sqrt(NPV.var))
  #based on objective Bayesian method
  B<-10000
  Se_seq<-rbeta(B,s1+1/2,n1-s1+1/2)
  FPR_seq<-rbeta(B,r1+1/2,n0-r1+1/2)
  p_seq<-rbeta(B,x+1/2,n-x+1/2)
  PPV_seq<-p_seq*Se_seq/(p_seq*Se_seq+(1-p_seq)*(FPR_seq))
  NPV_seq<-(1-p_seq)*(1-FPR_seq)/((1-FPR_seq)*(1-p_seq)+(1-Se_seq)*p_seq)
  PPV.interval3<-c(quantile(PPV_seq,probs=alpha/2),quantile(PPV_seq,probs=1-alpha/2))
  NPV.interval3<-c(quantile(NPV_seq,probs=alpha/2),quantile(NPV_seq,probs=1-alpha/2))
  output$Estimate_relevant_population<-c(PPV,NPV)
  output$PPV.interval1_relevant_population<-PPV.interval1
  output$PPV.interval2_relevant_population<-PPV.interval2
  output$PPV.interval3_relevant_population<-PPV.interval3
  output$NPV.interval1_relevant_population<-NPV.interval1
  output$NPV.interval2_relevant_population<-NPV.interval2
  output$NPV.interval3_relevant_population<-NPV.interval3
  return(output)
}
