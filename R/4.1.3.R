#' Compute Sensitivity, Specificity and Predictive Values with Clustered Binary-Scale Data
#' @param data a 3-dimensional array
#'             each cluster is represented by a 2*2 table(a matrix)
#'             the third dimension is the number of cluster(the number of matrices)
#' @param alpha significance level
#' @return estimation of Se,Sp,PPV,NPV,their variances and confidence intervals(wald's interval)
#' #example
#' #matrix1 <- matrix(c(22, 3, 2, 3), nrow = 2, byrow = TRUE)
#' #matrix2 <- matrix(c(56, 6, 23, 78), nrow = 2, byrow = TRUE)
#' #data <- array(c(matrix1, matrix2), dim = c(2, 2, 2))
#' #roc.clusteredbinary(data,0.05)
#' @export
#---------------------4.1.3 Sensitivity, Specificity and Predictive Values with Clustered Binary-Scale Data
roc.clusteredbinary<-function(data,alpha){
  result<-list(Se=NULL,Se.var=NULL,Sp=NULL,Sp.var=NULL,PPV=NULL,PPV.var=NULL,NPV=NULL,NPV.var=NULL,Se.interval=NULL,Sp.interval=NULL,PPV.interval=NULL,NPV.interval=NULL)
  I<-dim(data)[3]
  for(i in 1:I){
    if(dim(data)[1]!=2||dim(data)[2]!=2){
      stop("each cluster should represented by an 2*2 table")
    }
  }
  N_Se<-rep(NA,I)
  for(i in 1:I){
    N_Se[i]<-sum(data[1,,i])
  }
  N_Sp<-rep(NA,I)
  for(i in 1:I){
    N_Sp[i]<-sum(data[2,,i])
  }
  N_PPV<-rep(NA,I)
  for(i in 1:I){
    N_PPV[i]<-sum(data[,1,i])
  }
  N_NPV<-rep(NA,I)
  for(i in 1:I){
    N_NPV[i]<-sum(data[,2,i])
  }

  Se_seq<-data[1,1,]/(data[1,1,]+data[1,2,])
  Sp_seq<-data[2,2,]/(data[2,1,]+data[2,2,])
  PPV_seq<-data[1,1,]/(data[1,1,]+data[2,1,])
  NPV_seq<-data[2,2,]/(data[1,2,]+data[2,2,])
  Se<-sum(Se_seq*N_Se)/sum(N_Se)
  Sp<-sum(Sp_seq*N_Sp)/sum(N_Sp)
  PPV<-sum(PPV_seq*N_PPV)/sum(N_PPV)
  NPV<-sum(NPV_seq*N_NPV)/sum(N_NPV)
  Se.var<-sum((N_Se/mean(N_Se))^2*(Se_seq-Se)^2)/(I*(I-1))
  Sp.var<-sum((N_Sp/mean(N_Sp))^2*(Sp_seq-Sp)^2)/(I*(I-1))
  PPV.var<-sum((N_PPV/mean(N_PPV))^2*(PPV_seq-PPV)^2)/(I*(I-1))
  NPV.var<-sum((N_NPV/mean(N_NPV))^2*(NPV_seq-NPV)^2)/(I*(I-1))
  z<-qnorm(1-alpha/2)
  result$Se<-Se
  result$Se.var<-Se.var
  result$Sp<-Sp
  result$Sp.var<-Sp.var
  result$PPV<-PPV
  result$PPV.var<-PPV.var
  result$NPV<-NPV
  result$NPV.var<-NPV.var
  result$Se.interval<-c(Se-z*sqrt(Se.var),Se+z*sqrt(Se.var))
  result$Sp.interval<-c(Sp-z*sqrt(Sp.var),Sp+z*sqrt(Sp.var))
  result$PPV.interval<-c(PPV-z*sqrt(PPV.var),PPV+z*sqrt(PPV.var))
  result$NPV.interval<-c(NPV-z*sqrt(NPV.var),NPV+z*sqrt(NPV.var))
  return(result)
}
