#' Draw an empirical ROC curve and calculate AUC for ordinal data
#' Also compute variance and CI for AUC
#' @param data a matrix with two rows.
#'             The first row is test results for diseased subjects
#'             The second row is test results for nondiseased subjects
#' @param alpha significance level
#' @return estimation of AUC and its variance and confidence interval(using Hanley's and Delong's method)
#' #example
#' #data<-matrix(c(1,2,3,14,42,38,25,15,19,4),nrow=2,byrow=TRUE)
#' #roc.ordinal.emp(data,0.05)
#' @export
#-----4.2.1 Draw a scatter plot
roc.ordinal.emp<-function(data, alpha){

  result<-list( AUC=NULL, AUC_VAR=NULL, AUC_CI=NULL)

  if(dim(data)[1]!=2){
    stop("The row for table data should be 2")
  }

  K<-dim(data)[2]
  n1<-sum(data[1,])
  n0<-sum(data[2,])

  #---- draw empirical ROC curve ----------
  se<-c(rep(NA, length = K),0)
  FPR<-c(rep(NA, length = K),0)
  for(i in 1:K){
    se[i]= sum(data[1,i:K])/n1
    FPR[i] = sum(data[2,i:K])/n0
  }
  plot(FPR,se, type="b", xlab="FPR", ylab="Sensitivity", col="green",main=paste("Empirical ROC curve \nfor Ordinal Data with ", K, "Categories"))

  #--- area under empirical ROC curve ------
  tmpM<-matrix(NA, nrow=K, ncol=K)
  for(i in 1:K){
    for(j in 1:K){
      if(j > i){
        tmpM[i,j]<-0
      }
      if(i == j){
        tmpM[i,j]<-0.5
      }
      if(j < i){
        tmpM[i,j]<-1
      }
    }
  }
  tmpM2<-matrix(NA, ncol=K, nrow=K)
  for(i in 1:K){
    for(j in 1:K){
      tmpM2[i,j]<-data[1,i]*data[2,j]*tmpM[i,j]
    }
  }
  auc<-sum(tmpM2)/n0/n1

  #----variance of AUC Hanley & McNeil---------
  Q1<-auc/(2-auc)
  Q2<-2*auc^2/(1+auc)
  auc.var_HM<-(auc*(1-auc) + (n1-1)*(Q1-auc^2) + (n0-1)*(Q2-auc^2))/n0/n1

  #---variance of AUC by Delong, Delong & Clarke-Pearson--
  V1<-rep(NA, length=K)
  V0<-rep(NA, length=K)
  tmpM.V1<-matrix(NA, ncol=K, nrow=K)
  tmpM.V0<-matrix(NA, ncol=K, nrow=K)
  for(j in 1:K){
    tmpM.V1[,j]<-data[2,j]*tmpM[,j]
  }
  V1<-apply(tmpM.V1, 1, sum)/n0
  for(i in 1:K){
    tmpM.V0[i,]<-data[1,i]*tmpM[i,]
  }
  V0<-apply(tmpM.V0, 2, sum)/n1
  SS.T1<-sum(data[1,]*(V1-auc)^2)/(n1-1)
  SS.T0<-sum(data[2,]*(V0-auc)^2)/(n0-1)
  auc.var_DDC<-SS.T1/n1 + SS.T0/n0


  #----- calculate CI for AUC ---------------
  z<-qnorm(1-alpha/2)
  auc.ci_HM <-c(auc-z*sqrt(auc.var_HM), auc+z*sqrt(auc.var_HM))
  auc.ci_DDC<-c(auc-z*sqrt(auc.var_DDC), auc+z*sqrt(auc.var_DDC))

  result$AUC <- auc
  result$AUC_VAR<- rbind(auc.var_HM,auc.var_DDC)
  result$AUC_CI<- rbind(auc.ci_HM,auc.ci_DDC)
  return(result)
}
