#' using nonparametric method to estimate full area and partial area under ROC curve
#' also estimate the variance of full area
#' @param data0 test results of nondiseased subjects
#' @param data1 test results of diseased subjects
#' @param p the partial area between the FPRs e1 and e2, which correspond to the p-th and s-th ordinal test results
#' @param s the partial area between the FPRs e1 and e2, which correspond to the p-th and s-th ordinal test results
#' @return estimation of full area using MW estimate and Delong's estimate and their variances,
#' estimation of partial area
#' #example
#' #data0<-c(38,25,15,19,4)
#' #data1<-c(1,2,3,14,42)
#' #roc.nonparametric(data0,data1,2,4)
#' @export
roc.nonparametric<-function(data0,data1,p,s){
  generate<-function(data){
    result<-c()
    for(i in 1:length(data)){
      result<-c(result,rep(i,data[i]))
    }
    return(result)
  }
  x<-generate(data0)
  y<-generate(data1)
  #------MW area estimate------
  AMW<-0
  n0<-length(x)
  n1<-length(y)
  for(i in 1:n0){
    for (j in 1:n1){
      if(x[i]<y[j]) {AMW<-AMW+1}
      if(x[i]==y[j]) {AMW<-AMW+1/2}
    }
  }
  AMW<-AMW/(n0*n1)
  Q1<-AMW/(2-AMW)
  Q2<-2*AMW^2/(1+AMW)
  AMW.var<-(AMW*(1-AMW)+(n1-1)*(Q1-AMW^2)+(n0-1)*(Q2-AMW^2))/(n0*n1)

  #------Delong area estimate------
  a<-rep(0,times=length(x))
  b<-rep(0,times=length(y))
  for(j in 1:n0){
    for(i in 1:n1)
    {
      if(y[i]>x[j])
        a[j]<-a[j]+1/n1
      if(y[i]==x[j])
        a[j]<-a[j]+1/(2*n1)
    }
  }
  for(j in 1:n1){
    for(i in 1:n0)
    {
      if(y[j]>x[i])
        b[j]<-b[j]+1/n0
      if(y[j]==x[i])
        b[j]<-b[j]+1/(2*n0)
    }
  }
  V10<-mean(b)
  V01<-mean(a)
  ADL<-(V10+V01)/2
  S10<-0
  S01<-0
  for(i in 1:n1){
    S10<-S10+1/(n1-1)*(b[i]-ADL)^2
  }
  for(j in 1:n0){
    S01<-S01+1/(n0-1)*(a[j]-ADL)^2
  }
  ADL.var<-1/n1*S10+1/n0*S01
  partialarea<-0
  for(i in 1:n1){
    for(j in 1:n0){
      if((y[i]>x[j])&&(x[j]>=p)&&(x[j]<=s))
        partialarea<-partialarea+1/(n0*n1)
      if((y[i]==x[j])&&(x[j]>=p)&&(x[j]<=s))
        partialarea<-partialarea+1/(2*n0*n1)
    }
  }
  result<-list(AMW=NULL,ADL=NULL,AMW.var=NULL,ADL.var=NULL,partialarea=NULL)
  result$AMW<-AMW
  result$ADL<-ADL
  result$AMW.var<-AMW.var
  result$ADL.var<-ADL.var
  result$partialarea<-partialarea
  return(result)
}
