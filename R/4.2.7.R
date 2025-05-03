#' using nonparametric method to estimate AUC and its variance
#' @param data the test result.It is a matrix consisting of even rows
#' each cluster corresponds to two rows
#' The first row is the ordinal test results for diseased units in that cluster
#' The second row is the ordinal test results for undiseased units in that cluster
#' The cols represent ordinal test results.
#' @param alpha 1-alpha is the desired confidence level
#' @param B number of bootstrap replications
#' @return estimation of AUC and its variance
#' #example
#' #data<-matrix(c(0,0,0,3,1,2,0,1,0,0,0,0,1,0,0,1,1,1,0,0,0,0,2,0,1,0,1,0,1,0,0,0,2,0,3,0,3,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,2,0,2,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,2,0,0,1,0,0,0,2,0,0,0,1,0,0),nrow=20,byrow=TRUE)
#' #roc.clustered(data,0.05,1000)
#' @export
roc.clustered<-function(data,alpha,B){
  if((dim(data)[1]%%2)!=0){
    stop("Each cluster should contain 2 rows")
  }
  psi<-function(x,y){
    ifelse(x>y,1,ifelse(x<y,0,0.5))
  }
  result<-list(Ac=NULL,Ac.var=NULL,CI=NULL)
  K<-dim(data)[2]
  I<-dim(data)[1]/2
  odd<-seq(1,dim(data)[1],by=2)
  yes_matrix<-data[odd,]
  even<-seq(2,dim(data)[1],by=2)
  no_matrix<-data[even,]
  Ac<-0
  for(i in 1:I){
    for(i_ in 1:I){
      if(sum(yes_matrix[i,])>0&&sum(no_matrix[i_,])>0){
        tmp_yes<-rep(seq_along(yes_matrix[i,]),times=yes_matrix[i,])
        tmp_no<-rep(seq_along(no_matrix[i_,]),times=no_matrix[i_,])
        for(j in 1:length(tmp_yes)){
          for(k in 1:length(tmp_no)){
            Ac<-Ac+psi(tmp_yes[j],tmp_no[k])
          }
        }
      }
    }
  }
  Ac<-Ac/(sum(yes_matrix)*sum(no_matrix))
  #对固定i，定义V
  I10<-sum(rowSums(yes_matrix!=0)>0)
  I01<-sum(rowSums(no_matrix!=0)>0)
  n1<-rowSums(yes_matrix)
  n0<-rowSums(no_matrix)
  S10<-0
  S01<-0
  S11<-0
  #直接计算各项之和
  V10_sum<-rep(0,I)
  for(i in 1:I){
    for(i_ in 1:I){
      if(sum(yes_matrix[i,])>0&&sum(no_matrix[i_,])>0){
        tmp_yes<-rep(seq_along(yes_matrix[i,]),times=yes_matrix[i,])
        tmp_no<-rep(seq_along(no_matrix[i_,]),times=no_matrix[i_,])
        for(j in 1:length(tmp_yes)){
          for(k in 1:length(tmp_no)){
            V10_sum[i]<-V10_sum[i]+psi(tmp_yes[j],tmp_no[k])
          }
        }
      }
    }
  }
  V10_sum<-V10_sum/sum(n0)
  #直接计算各项之和
  V01_sum<-rep(0,I)
  for(i in 1:I){
    for(i_ in 1:I){
      if(sum(yes_matrix[i,])>0&&sum(no_matrix[i_,])>0){
        tmp_yes<-rep(seq_along(yes_matrix[i,]),times=yes_matrix[i,])
        tmp_no<-rep(seq_along(no_matrix[i_,]),times=no_matrix[i_,])
        for(j in 1:length(tmp_yes)){
          for(k in 1:length(tmp_no)){
            V01_sum[i]<-V01_sum[i]+psi(tmp_yes[j],tmp_no[k])
          }
        }
      }
    }
  }
  V01_sum<-V01_sum/sum(n1)
  S10<-sum((V10_sum-n1*Ac)^2)*I10/(I10-1)/sum(n1)
  S01<-sum((V01_sum-n0*Ac)^2)*I01/(I01-1)/sum(n0)
  S11<-sum((V10_sum-n1*Ac)*(V01_sum-n0*Ac))*I/(I-1)
  Ac.var<-S10/sum(n1)+S01/sum(n0)+S11/sum(n0)/sum(n1)
  result$Ac<-Ac
  result$Ac.var<-Ac.var
  #calculate confidence interval based on Bootstrap percentile method
  Ac_boot<-rep(0,B)
  for(iter in 1:B){
    indices<-sample(1:I,replace=TRUE)
    yes_matrix<-matrix(NA,nrow=I,ncol=K)
    no_matrix<-matrix(NA,nrow=I,ncol=K)
    for(i in 1:I){
      yes_matrix[i,]<-data[2*indices[i]-1,]
      no_matrix[i,]<-data[2*indices[i],]
    }
    Ac_boot[iter]<-0
    for(i in 1:I){
      for(i_ in 1:I){
        if(sum(yes_matrix[i,])>0&&sum(no_matrix[i_,])>0){
          tmp_yes<-rep(seq_along(yes_matrix[i,]),times=yes_matrix[i,])
          tmp_no<-rep(seq_along(no_matrix[i_,]),times=no_matrix[i_,])
          for(j in 1:length(tmp_yes)){
            for(k in 1:length(tmp_no)){
              Ac_boot[iter]<-Ac_boot[iter]+psi(tmp_yes[j],tmp_no[k])
            }
          }
        }
      }
    }
    Ac_boot[iter]<-Ac_boot[iter]/(sum(yes_matrix)*sum(no_matrix))
  }
  result$CI<-c(quantile(Ac_boot,probs=alpha/2),quantile(Ac_boot,probs=1-alpha/2))
  return(result)
}

