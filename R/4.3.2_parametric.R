#' using a Box-Cox transformation to transform continuous data into binormal data
#' @param data0 continuous test results of nondiseased individuals
#' @param data1 continuous test results of diseased individuals
#' @return the optimal parameter value for the Box-Cox transformation,the binormal parameters after transformation and corresponding AUC
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #BoxCox(data0,data1)
#' @export
BoxCox<-function(data0,data1){
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
  result<-nlminb(start=1,loglikelihood)
  lambda<-result$par
  if(abs(lambda)>=1e-4){
    data0_transformed<-(data0^lambda-1)/lambda
    data1_transformed<-(data1^lambda-1)/lambda
  }
  else{
    data0_transformed<-log(data0)
    data1_transformed<-log(data1)
  }
  a<-(mean(data1_transformed)-mean(data0_transformed))/sqrt(var(data1_transformed)*(n1-1)/n1)
  b<-sqrt(var(data0_transformed)*(n0-1)/n0)/sqrt(var(data1_transformed)*(n1-1)/n1)
  area.full<-pnorm(a/sqrt(1+b^2))
  # 输出结果
  output<-list(lambda=NULL,a=NULL,b=NULL,area.full=NULL)
  output$lambda<-lambda
  output$a<-a
  output$b<-b
  output$area.full<-area.full
  return(output)
}

#' Function to calculating binormal parameters and plot the smooth ROC curve for continuous data
#' @param data0 continuous test results of nondiseased individuals
#' @param data1 continuous test results of diseased individuals
#' @return the optimal parameter value for the Box-Cox transformation,binormal parameters and their estimated covariance matrix,and full area
#' #example
#' #data1<-c(140,1087,230,183,1256,700,16,800,253,740,126,153,283,90,303,193,76,1370,543,913,230,463,60,509,576,671,80,490,156,356,350,323,1560,120,216,443,523,76,303,353,206)
#' #data0<-c(136,286,281,23,200,146,220,96,100,60,17,27,126,100,253,70,40,6,46)
#' #roc.continuous.smooth(data0,data1)
#' @export
roc.continuous.smooth<-function(data0,data1){
  if(sum(data0<=0)+sum(data1<=0)>0){
    stop("The test results should be positive numbers")
  }
  output<-list(lambda=NULL,a=NULL,b=NULL,cov_matrix=matrix(rep(NA,4),nrow=2),area.full=NULL)
  lambda<-BoxCox(data0,data1)$lambda
  a<-BoxCox(data0,data1)$a
  b<-BoxCox(data0,data1)$b
  n0<-length(data0)
  n1<-length(data1)
  if(n0>=20&&n1>=20){
    output$lambda<-lambda
    output$a<-a
    output$b<-b
    output$cov_matrix[1,1]<-(n0*(a^2+2)+2*n1*b^2)/(2*n0*n1)
    output$cov_matrix[2,2]<-(n1+n0)*b^2/(2*n0*n1)
    output$cov_matrix[1,2]<-a*b/(2*n1)
    output$cov_matrix[2,1]<-a*b/(2*n1)
    output$area.full<-BoxCox(data0,data1)$area.full
  }
  #样本量较小时用Jackknife刀切法估计方差
  else{
   a<-BoxCox(data0,data1)$a
   b<-BoxCox(data0,data1)$b
   n0<-length(data0)
   n1<-length(data1)
   a_seq<-rep(NA,n0+n1)
   b_seq<-rep(NA,n0+n1)
   for(i in 1:n0){
     data0_new<-data0[-i]
     a_seq[i]<-BoxCox(data0_new,data1)$a
     b_seq[i]<-BoxCox(data0_new,data1)$b
   }
   for(j in 1:n1){
    data1_new<-data1[-j]
    a_seq[j+n0]<-BoxCox(data0,data1_new)$a
    b_seq[j+n0]<-BoxCox(data0,data1_new)$b
    }
    output$lambda<-lambda
    output$a<-a
    output$b<-b
    output$cov_matrix[1,1]<-(n0+n1-1)*sum((a_seq-mean(a_seq))^2)/(n0+n1)
    output$cov_matrix[2,2]<-(n0+n1-1)*sum((b_seq-mean(b_seq))^2)/(n0+n1)
    output$cov_matrix[1,2]<-(n0+n1-1)*sum((a_seq-mean(a_seq))*(b_seq-mean(b_seq)))/(n0+n1)
    output$cov_matrix[2,1]<-(n0+n1-1)*sum((a_seq-mean(a_seq))*(b_seq-mean(b_seq)))/(n0+n1)
    output$area.full<-BoxCox(data0,data1)$area.full
  }
  #绘制曲线
  x<-seq(0,1,length=10000)
  y<-rep(0,length=10000)
  for(i in 1:10000){
    y[i]<-pnorm(a+b*qnorm(x[i]))
  }
  plot(x,y,type="l",main="Estimation of Smooth ROC Curve",xlab="FPR",ylab="TPR",col="purple")
  return(output)
}





