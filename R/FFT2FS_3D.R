FFT2FS_3D <-
function(A){
  ##calculate the dimension of the array
  n1 <- dim(A)[1]  
  n2 <- dim(A)[2]
  n3 <- dim(A)[3]
  
  ##calculate the FFT of the 3-D array
  FFT_A <- fft(A)/(n1*n2*n3)
  
  FFT2FS_0 <- function(X) {    
    ## partition the X matrix into four parts
    F11 <- X[1,1]
    F12 <- X[1,-1]
    F21 <- X[-1,1]
    F22 <- X[-1,-1]
    
    C11 <- as.matrix(Re(F11),nrow=1)
    ## Calculation of C12
    C12_r <- matrix(Re(F12),nrow=1)
    C12_i <- matrix(-Im(F12),nrow=1)
    colnames(C12_r) <- seq(1,2*(n2-1),by=2)
    colnames(C12_i) <- seq(2,2*(n2-1),by=2)
    C12 <- cbind(C12_r,C12_i)   
    C12 <- C12[1, as.character(sort(as.numeric(colnames(C12))))]
    ## Calculation of C21
    C21_r <- matrix(Re(F21),ncol=1)
    C21_i <- matrix(-Im(F21),ncol=1)
    rownames(C21_r) <- seq(1,2*(n1-1),by=2)
    rownames(C21_i) <- seq(2,2*(n1-1),by=2)
    C21 <- rbind(C21_r,C21_i)
    C21 <- C21[as.character(sort(as.numeric(rownames(C21)))),1]
    ## Calculation of C22
    C22_1 <- Re(F22)
    C22_2 <- -Im(F22)
    C22_3 <- -Im(F22)
    C22_4 <- -Re(F22)
    # Calculation of C22_odd
    colnames(C22_1) <- seq(1,2*(n2-1),by=2)
    colnames(C22_2) <- seq(2,2*(n2-1),by=2)
    C22_odd <- cbind(C22_1,C22_2)
    C22_odd <- C22_odd[,as.character(sort(as.numeric(colnames(C22_odd))))]
    rownames(C22_odd) <- seq(1,2*(n1-1),by=2)
    
    # Calculate of the C22_even
    colnames(C22_3) <- seq(1,2*(n2-1),by=2)
    colnames(C22_4) <- seq(2,2*(n2-1),by=2)
    C22_even <- cbind(C22_3,C22_4)
    C22_even <- C22_even[,as.character(sort(as.numeric(colnames(C22_even))))]
    rownames(C22_even) <- seq(2,2*(n1-1),by=2)
    
    #Build the C22 matrix
    C22 <- rbind(C22_odd,C22_even)
    C22 <- C22[as.character(sort(as.numeric(rownames(C22)))),]
    
    ###calculations of the C matrix
    C_1 <- matrix(c(C11,C12),nrow=1)
    C_2 <- cbind(C21,C22)
    C <- rbind(C_1,C_2)
    ### return the C matrix
    return(C)
  }
  
  
  ###FFT to FS for the cos frequency
  FFT2FS_cos <- function(X) {
    
    ## partition the X matrix into four parts
    F11 <- X[1,1]
    F12 <- X[1,-1]
    F21 <- X[-1,1]
    F22 <- X[-1,-1]
    
    C11 <- as.matrix(Re(F11),nrow=1)
    ## Calculation of C12
    C12_r <- matrix(Re(F12),nrow=1)
    C12_i <- matrix(-Im(F12),nrow=1)
    colnames(C12_r) <- seq(1,2*(n2-1),by=2)
    colnames(C12_i) <- seq(2,2*(n2-1),by=2)
    C12 <- cbind(C12_r,C12_i)   
    C12 <- C12[1, as.character(sort(as.numeric(colnames(C12))))]
    ## Calculation of C21
    C21_r <- matrix(Re(F21),ncol=1)
    C21_i <- matrix(-Im(F21),ncol=1)
    rownames(C21_r) <- seq(1,2*(n1-1),by=2)
    rownames(C21_i) <- seq(2,2*(n1-1),by=2)
    C21 <- rbind(C21_r,C21_i)
    C21 <- C21[as.character(sort(as.numeric(rownames(C21)))),1]
    ## Calculation of C22
    C22_1 <- Re(F22)
    C22_2 <- -Im(F22)
    C22_3 <- -Im(F22)
    C22_4 <- -Re(F22)
    # Calculation of C22_odd
    colnames(C22_1) <- seq(1,2*(n2-1),by=2)
    colnames(C22_2) <- seq(2,2*(n2-1),by=2)
    C22_odd <- cbind(C22_1,C22_2)
    C22_odd <- C22_odd[,as.character(sort(as.numeric(colnames(C22_odd))))]
    rownames(C22_odd) <- seq(1,2*(n1-1),by=2)
    
    # Calculate of the C22_even
    colnames(C22_3) <- seq(1,2*(n2-1),by=2)
    colnames(C22_4) <- seq(2,2*(n2-1),by=2)
    C22_even <- cbind(C22_3,C22_4)
    C22_even <- C22_even[,as.character(sort(as.numeric(colnames(C22_even))))]
    rownames(C22_even) <- seq(2,2*(n1-1),by=2)
    
    #Build the C22 matrix
    C22 <- rbind(C22_odd,C22_even)
    C22 <- C22[as.character(sort(as.numeric(rownames(C22)))),]
    
    ###calculations of the C matrix
    C_1 <- matrix(c(C11,C12),nrow=1)
    C_2 <- cbind(C21,C22)
    C <- rbind(C_1,C_2)
    ### return the C matrix
    return(C)
  }
  
  ###FFT to FS for the sin frequency
  FFT2FS_sin <- function(X) {
    ##calculate the dimension of the matrix
    
    ## partition the X matrix into four parts
    F11 <- X[1,1]
    F12 <- X[1,-1]
    F21 <- X[-1,1]
    F22 <- X[-1,-1]
    
    C11 <- as.matrix(-Im(F11),nrow=1)
    ## Calculation of C12
    C12_r <- matrix(-Re(F12),nrow=1)
    C12_i <- matrix(-Im(F12),nrow=1)
    colnames(C12_i) <- seq(1,2*(n2-1),by=2)
    colnames(C12_r) <- seq(2,2*(n2-1),by=2)
    C12 <- cbind(C12_i,C12_r)   
    C12 <- C12[1, as.character(sort(as.numeric(colnames(C12))))]
    
    ## Calculation of C21
    C21_r <- matrix(-Re(F21),ncol=1)
    C21_i <- matrix(-Im(F21),ncol=1)
    rownames(C21_i) <- seq(1,2*(n1-1),by=2)
    rownames(C21_r) <- seq(2,2*(n1-1),by=2)
    C21 <- rbind(C21_i,C21_r)
    C21 <- C21[as.character(sort(as.numeric(rownames(C21)))),1]
    
    ## Calculation of C22
    C22_1 <- -Im(F22)
    C22_2 <- -Re(F22)
    C22_3 <- -Re(F22)
    C22_4 <-  Im(F22)
    
    # Calculation of C22_odd
    colnames(C22_1) <- seq(1,2*(n2-1),by=2)
    colnames(C22_2) <- seq(2,2*(n2-1),by=2)
    C22_odd <- cbind(C22_1,C22_2)
    C22_odd <- C22_odd[,as.character(sort(as.numeric(colnames(C22_odd))))]
    rownames(C22_odd) <- seq(1,2*(n1-1),by=2)
    
    # Calculate of the C22_even
    colnames(C22_3) <- seq(1,2*(n2-1),by=2)
    colnames(C22_4) <- seq(2,2*(n2-1),by=2)
    C22_even <- cbind(C22_3,C22_4)
    C22_even <- C22_even[,as.character(sort(as.numeric(colnames(C22_even))))]
    rownames(C22_even) <- seq(2,2*(n1-1),by=2)
    
    #Build the C22 matrix
    C22 <- rbind(C22_odd,C22_even)
    C22 <- C22[as.character(sort(as.numeric(rownames(C22)))),]
    
    ###calculations of the C matrix
    C_1 <- matrix(c(C11,C12),nrow=1)
    C_2 <- cbind(C21,C22)
    C <- rbind(C_1,C_2)
    ### return the C matrix
    return(C)
  }
  
  
  
  
  C_matrix <- array(0,dim=c(2*n1-1,2*n2-1,2*n3-1))
  ##mapping the coefficients at the w=0 plane
  C_matrix[,,1] <- FFT2FS_0(FFT_A[,,1])
  ## mapping the rest of the coefficients by using cos and sin frequencies
  for (i in 1:(n3-1)){
    C_matrix[,,2*i] <- FFT2FS_cos(FFT_A[,,i+1])
    C_matrix[,,2*i+1] <- FFT2FS_sin(FFT_A[,,i+1])
  }
  return(C_matrix)
  
}
