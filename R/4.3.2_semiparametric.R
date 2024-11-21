#' Function used to fit a smooth ROC curve.It is not important.
#' @param Y_star not important
#' @param k_star not important
#' @param l_star not important
#' @export
star<-function(Y_star,k_star,l_star){
  result<-rep(0,length=length(Y_star))
  for(r in 1:length(Y_star)){
    if(k_star[r]>0&&l_star[r]>0){
      result[r]<-2
    }
    else if(k_star[r]>0&&l_star[r]==0){
      result[r]<-1
    }
    else{
      result[r]<-0
    }
  }
  return(result)
}

#' Function used to fitting a smooth ROC curve.It is not important.
#' @param data0 continuous test result of undiseased subjects
#' @param data1 continuous test result of diseased subjects
#' @export
getconstant<-function(data0,data1){
  Y_star<-sort(unique(c(data0,data1)))
  k_star<-rep(0,length=length(Y_star))
  l_star<-rep(0,length=length(Y_star))
  for(r in 1:length(Y_star)){
    k_star[r]<-sum(data0==Y_star[r])
    l_star[r]<-sum(data1==Y_star[r])
  }
  Y<-c()
  for(r in 1:(length(Y_star)-2)){
    if(star(Y_star,k_star,l_star)[r]!=star(Y_star,k_star,l_star)[r+1]){
      Y<-c(Y,Y_star[r])
    }
    else if(star(Y_star,k_star,l_star)[r]==2||star(Y_star,k_star,l_star)[r+1]==2){
      Y<-c(Y,Y_star[r])
    }
  }
  k<-rep(0,length=length(Y)+1)
  l<-rep(0,length=length(Y)+1)
  for(r in 2:length(Y)){
    k[r]<-sum(data0>Y[r-1]&data0<=Y[r])
    l[r]<-sum(data1>Y[r-1]&data1<=Y[r])
  }
  k[1]<-sum(data0<=Y[1])
  l[1]<-sum(data1<=Y[1])
  k[length(Y)+1]<-sum(data0>Y[length(Y)])
  l[length(Y)+1]<-sum(data1>Y[length(Y)])
  output<-list(k=NULL,l=NULL)
  output$k<-k
  output$l<-l
  return(output)
}

#' Fitting a smooth ROC curve for continuous data using semiparametric method
#' @param data0 continuous test result for undiseased subjects
#' @param data1 continuous test result for diseased subjects
#' @return a smooth ROC curve
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.semiparametric(data0,data1)
#' @export
roc.semiparametric<-function(data0,data1){
  k<-getconstant(data0,data1)$k
  l<-getconstant(data0,data1)$l
  loglikelihood<-function(params){
    C<-params[1:(length(k)-1)]
    alpha0<-params[length(k)]
    alpha1<-params[length(k)+1]
    if(any(diff(C)<=0)||alpha1<=0){
      return(Inf)
    }
    result<-k[1]*log(pnorm(C[1]))+l[1]*log(pnorm(-alpha0+alpha1*C[1]))
    result<-result+k[length(k)]*log(1-pnorm(C[length(k)-1]))+l[length(l)]*log(1-pnorm(-alpha0+alpha1*C[length(l)-1]))
    for(r in 2:(length(k)-1)){
      result<-result+k[r]*log(pnorm(C[r])-pnorm(C[r-1]))+l[r]*log(pnorm(-alpha0+alpha1*C[r])-pnorm(-alpha0+alpha1*C[r-1]))
    }
    return(-result)
  }
  initial_values<-c(seq(-1,1,length=length(k)-1),1,1)
  result<-nlminb(start=initial_values,loglikelihood,control=list(iter.max=2000,eval.max=100000))
  print(result)
  estimation<-result$par
  alpha0<-estimation[length(estimation)-1]
  alpha1<-estimation[length(estimation)]
  x<-seq(0,1,by=0.01)
  y<-pnorm(alpha0+alpha1*qnorm(x))
  plot(x,y,xlab="FPR",ylab="TPR",type="l",main="Estimation of Smooth ROC Curve Using Semiparametric Method",col="brown")
}


