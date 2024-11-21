#' testing the hypothesis that the roc curve area or partial area is a specific value
#' @param A0 area or partial area under the null hypothesis
#' @param A the estimated value of area or partial area
#' @param variance estimated variance of area or partial area
#' @param alpha significance level
#' @return whether the null hypothesis is rejected
#' #example
#' #roc.testing(0.5,0.9297,0.00045369,0.05)
#' @export
roc.testing<-function(A0,A,variance,alpha){
  statistic<-(A-A0)/sqrt(variance)
  print(c("The test statistic is",statistic))
  if(abs(statistic)>qnorm(1-alpha/2)){
    print("The null hypothesis is rejected")
  }
  else{
    print("The null hypothesis is accepted")
  }
}
