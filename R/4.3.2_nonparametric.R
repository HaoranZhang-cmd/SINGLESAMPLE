#' function to transform the original data to normal data using Box-Cox transformation
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @return the transformed data and the optimal parameter value for the Box-Cox transformation
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.transformed(data0,data1)
#' @export
roc.transformed<-function(data0,data1){
  n0<-length(data0)
  n1<-length(data1)
  loglikelihood<-function(lambda){
    if(abs(lambda)>=1e-4){
      data0_transformed<-(data0^lambda-1)/lambda
      data1_transformed<-(data1^lambda-1)/lambda
    }
    else{
      data0_transformed<-log(data0)
      data1_transformed<-log(data1)
    }
    mean0<-mean(data0_transformed)
    var0<-var(data0_transformed)*(n0-1)/n0
    result0<-log(sqrt(var0))*(-n0)+sum(log(data0^(lambda-1))-(data0_transformed-mean0)^2/(2*var0))
    mean1<-mean(data1_transformed)
    var1<-var(data1_transformed)*(n1-1)/length(data1)
    result1<-log(sqrt(var1))*(-n1)+sum(log(data1^(lambda-1))-(data1_transformed-mean1)^2/(2*var1))
    return(-result0-result1)
  }
  lambda<-nlminb(start=1,loglikelihood)$par
  if(abs(lambda)>=1e-4){
    data0_transformed<-(data0^lambda-1)/lambda
    data1_transformed<-(data1^lambda-1)/lambda
  }
  else{
    data0_transformed<-log(data0)
    data1_transformed<-log(data1)
  }
  output<-list(data0_transformed=NULL,data1_transformed=NULL,lambda=NULL)
  output$data0_transformed<-data0_transformed
  output$data1_transformed<-data1_transformed
  output$lambda<-lambda
  return(output)
}

#' function to draw smooth ROC curve using biweight kernel
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @return a smooth ROC curve
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.biweight(data0,data1)
#' @export

roc.biweight<-function(data0,data1){
  new_data0<-roc.transformed(data0,data1)$data0_transformed
  new_data1<-roc.transformed(data0,data1)$data1_transformed
  int<-function(c,test,h){
    if(c>test+h){
      return(0)
    }
    else{
      begin<-max(c,test-h)
      end<-test+h
      return((end-test)^5/(5*h^4)-2*(end-test)^3/(3*h^2)-(begin-test)^5/(5*h^4)+2*(begin-test)^3/(3*h^2)+end-begin)
    }
  }
  n0<-length(new_data0)
  n1<-length(new_data1)
  h0<-0.9*min(sd(new_data0),(sort(new_data0)[ceiling(n0*0.75)]-sort(new_data0)[ceiling(n0*0.25)])/1.34)/(n0)^(1/5)
  h1<-0.9*min(sd(new_data1),(sort(new_data1)[ceiling(n1*0.75)]-sort(new_data1)[ceiling(n1*0.25)])/1.34)/(n1)^(1/5)
  F0<-function(c){
    result<-0
    for(j in 1:n0){
      result<-result+int(c,new_data0[j],h0)/(n0*h0)*15/16
    }
    return(result)
  }
  F1<-function(c){
    result<-0
    for(j in 1:n1){
      result<-result+int(c,new_data1[j],h1)/(n1*h1)*15/16
    }
    return(result)
  }
  cmin<-min(c(new_data0-h0,new_data1-h1))
  cmax<-max(c(new_data0+h0,new_data1+h1))
  c_seq<-seq(cmin,cmax,length=1000)
  x<-rep(0,1000)
  y<-rep(0,1000)
  for(i in 1:1000){
    x[i]<-F0(c_seq[i])
    y[i]<-F1(c_seq[i])
  }
  plot(x,y,type="l",col="green",xlab="FPR",ylab="TPR",main="Estimation of Smooth ROC Curve Using Biweight Kernel")
}


#' function to draw smooth ROC curve using gaussian kernel
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @return a smooth ROC curve
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.gaussian(data0,data1)
#' @export
roc.gaussian<-function(data0,data1){
  new_data0<-roc.transformed(data0,data1)$data0_transformed
  new_data1<-roc.transformed(data0,data1)$data1_transformed
  n0<-length(new_data0)
  n1<-length(new_data1)
  h0<-0.9*min(sd(new_data0),(sort(new_data0)[ceiling(n0*0.75)]-sort(new_data0)[ceiling(n0*0.25)])/1.34)/(n0)^(1/5)
  h1<-0.9*min(sd(new_data1),(sort(new_data1)[ceiling(n1*0.75)]-sort(new_data1)[ceiling(n1*0.25)])/1.34)/(n1)^(1/5)
  F0<-function(c){
    result<-0
    for(j in 1:n0){
      result<-result+pnorm((new_data0[j]-c)/h0)/n0
    }
    return(result)
  }
  F1<-function(c){
    result<-0
    for(j in 1:n1){
      result<-result+pnorm((new_data1[j]-c)/h1)/n1
    }
    return(result)
  }
  cmin<-min(c(new_data0-5*h0,new_data1-5*h1))
  cmax<-max(c(new_data0+5*h0,new_data1+5*h1))
  c_seq<-seq(cmin,cmax,length=1000)
  x<-rep(0,1000)
  y<-rep(0,1000)
  for(i in 1:1000){
    x[i]<-F0(c_seq[i])
    y[i]<-F1(c_seq[i])
  }
  plot(x,y,type="l",col="blue",xlab="FPR",ylab="TPR",main="Estimation of Smooth ROC Curve Using Gaussian Kernel")
}

#  using the method proposed by Zhou and Harezlar(2002)
#' kernel function proposed by Zhou and Herezlar
#' @param t calculate the kernel
#' @export
Epanechnikov_kernel<-function(t){
  if(t<(-1)){
    return(0)
  }
  else if(t>1){
    return(1)
  }
  return(3*t/4-t^3/4+1/2)
}

#' function to draw smooth ROC curve using gaussian kernel
#' @param data0 continuous test results of nondiseased subjects
#' @param data1 continuous test results of diseased subjects
#' @return a smooth ROC curve
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc_Zhou_Harezlar(data0,data1)
#' @export
roc_Zhou_Harezlar<-function(data0,data1){
  bandwidth<-matrix(c(0,0,0,0),nrow=2)
  n0<-length(data0)
  n1<-length(data1)
  new_data0<-roc.transformed(data0,data1)$data0
  new_data1<-roc.transformed(data0,data1)$data1
  bandwidth[1,1]<-0.9*min(sd(new_data0),(sort(new_data0)[ceiling(n0*0.75)]-sort(new_data0)[ceiling(n0*0.25)])/1.34)/(n0)^(1/5)
  bandwidth[1,2]<-(180*sqrt(pi)/(7*n0))^(1/3)*min(sd(data0),(sort(data0)[ceiling(length(data0)*0.75)]-sort(data0)[ceiling(length(data0)*0.25)])/1.34)
  bandwidth[2,1]<-0.9*min(sd(new_data1),(sort(new_data1)[ceiling(n1*0.75)]-sort(new_data1)[ceiling(n1*0.25)])/1.34)/(n1)^(1/5)
  bandwidth[2,2]<-(180*sqrt(pi)/(7*n1))^(1/3)*min(sd(data1),(sort(data1)[ceiling(length(data1)*0.75)]-sort(data1)[ceiling(length(data1)*0.25)])/1.34)
  #绘制曲线
  #第一条
  cut1<-seq(-max(c(new_data0,new_data1)),max(c(new_data0,new_data1)),length=1000)
  x1<-rep(0,length=length(cut1))
  y1<-rep(0,length=length(cut1))
  for(i in 1:length(x1)){
    for(j in 1:n0){
      x1[i]<-x1[i]+Epanechnikov_kernel((cut1[i]-new_data0[j])/bandwidth[1,1])/n0
    }
    for(j in 1:n1){
      y1[i]<-y1[i]+Epanechnikov_kernel((cut1[i]-new_data1[j])/bandwidth[2,1])/n1
    }
  }
  plot(y1,x1,xlab="FPR",ylab="TPR",xlim=c(0,1),ylim=c(0,1),type="l",col="purple",main="Estimation of Smooth ROC Curve Using Zou's Bandwidth")
  #第二条
  cut2<-seq(-max(c(data0,data1)),max(c(data0,data1)),length=1000)
  x2<-rep(0,length=length(cut2))
  y2<-rep(0,length=length(cut2))
  for(i in 1:length(x2)){
    for(j in 1:n0){
      x2[i]<-x2[i]+Epanechnikov_kernel((cut2[i]-data0[j])/bandwidth[1,2])/n0
    }
    for(j in 1:n1){
      y2[i]<-y2[i]+Epanechnikov_kernel((cut2[i]-data1[j])/bandwidth[2,2])/n1
    }
  }
  plot(y2,x2,xlab="FPR",ylab="TPR",xlim=c(0,1),ylim=c(0,1),type="l",col="pink",main="Estimation of Smooth ROC Curve Using Bowman's Bandwidth")
}
