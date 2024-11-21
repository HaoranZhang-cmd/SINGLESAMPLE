#' draw an empirical ROC curve for continuous data
#' @param data0 test results of undiseased subjects
#' @param data1 test results of diseased subjects
#' @return an empirical ROC curve
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.continuous.emp(data0,data1)
#' @export
roc.continuous.emp<-function(data0,data1){
  n0<-length(data0)
  n1<-length(data1)
  divide<-sort(c(data0,data1))
  c<-seq(divide[1]-1,divide[length(divide)]+1,by=0.1)
  x<-rep(0,length=length(c))
  y<-rep(0,length=length(c))
  for(i in 1:length(c)){
    for(j in 1:n0){
      x[i] <- x[i] + (data0[j] > c[i])/n0
    }
    for(j in 1:n1){
      y[i] <- y[i] + (data1[j] > c[i])/n1
    }
  }
  plot(x,y,type="l",xlab="FPR",ylab="Se",xlim=c(0,1),ylim=c(0,1),main="Empirical ROC Curve",col="green")
}
