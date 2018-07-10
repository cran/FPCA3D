FPCA_3D_score <-
function(X,prop){
  ##read the dimension of the array
  m <- dim(X)[1]
  n <- dim(X)[2]
  t <- dim(X)[3]
  k <- dim(X)[4]
  ####calculate the fourier series for each image and combine then into a big matrix
  FC <- matrix(rep(0,(2*m-1)*(2*n-1)*(2*t-1)*k),ncol=k)
  for (i in 1:k){
    FC[,i] <-  as.vector(FFT2FS_3D(X[,,,i]))
  }
  
  
  ####calculate the svd of the C matrix
  C_svd <- svd(FC/sqrt(k));
  C_value <- (C_svd$d)^2;
  C_prop = cumsum(C_value)/sum(C_value);
  C_ind = which(C_prop > prop)[1];
  C_vector <- C_svd$u;
  C_score <- t(FC) %*% C_vector;
  colnames(C_score) = sapply(1:length(C_value),function(x){paste("FPC_",x,sep="")});
  rlt = as.matrix(C_score[,c(1:C_ind)]);
  return(rlt)
}
