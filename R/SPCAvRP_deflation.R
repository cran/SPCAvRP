SPCAvRP_deflation <- function( data               # either the data matrix or the sample covariance matrix
                             , cov = FALSE        # TRUE if data is given as a sample covariance matrix 
                             , s                  # the dimension of the eigenspace (the number of principal components)
                             , l                  # sparsity level of the eigenvectors (array of lenght s) 
                             , d = 10             # dimension of the random projections (if k is known, use d = k)
                             , A = 300            # number of projections to aggregate 
                             , B = 100            # number of projections in a group
                             , center_data = TRUE # TRUE if the data matrix should be centered 
)
  # output : a list of two elements
  #             output$vector : a matrix of dimension pxs with its columns as estimated eigenvectors of sparsity level l 
  #             output$value  : an array of dimesnin s with estimated eigenvalues
{
  if( !cov & center_data ){
    data <- scale(data, center_data, FALSE)
  }
  p <- ncol(data)
  p_cond <- (p <= 2*d*sqrt(A*B))
  if( p_cond & cov == FALSE){
    data <- 1/nrow(data)*crossprod(data)
    cov <- TRUE
  }else if( !p_cond & cov == FALSE){
    print('SPCAvRP: since p > 2*d*sqrt(A*B), computation of p x p sample covariance is avoided')
  }
  
  lambda_hat <- rep(0,s)
  v_hat <- matrix(rep(0,s*p), nrow = p)
  v_hat_supp <- NULL
  spca_1 <- SPCAvRP(data, cov, l[1], d, A, B, center_data = FALSE) # can be any SPCA algorithm
  v_hat[,1] <- spca_1$vector
  lambda_hat[1] <- spca_1$value
  if(s>1){
    v_hat_supp[[1]] <-  which(v_hat[,1]!=0, arr.ind = T)
    Hdata <- data
    for(i in 2:s){
      H <- diag(length(v_hat_supp[[i-1]])) - tcrossprod(v_hat[v_hat_supp[[i-1]],i-1])
      if(cov == TRUE){
        Hdata[v_hat_supp[[i-1]],v_hat_supp[[i-1]]] <- H%*%Hdata[v_hat_supp[[i-1]],v_hat_supp[[i-1]]]%*%H
      }else{
        Hdata[,v_hat_supp[[i-1]]] <- Hdata[,v_hat_supp[[i-1]]]%*%H
      }
      spca_i <- SPCAvRP(Hdata, cov, l[i], d, A, B, center_data = FALSE) # can be any SPCA algorithm
      v_hat_supp[[i]] <- which(spca_i$vector!=0, arr.ind = T) 
      V <- NULL
      for(j in seq(i-1,1,-1)){
        u <- v_hat[v_hat_supp[[i]],j]
        if(norm(u,"2")!=0){V <- cbind(V,u/norm(u,"2"))}
      }
      if(is.null(V)){
        Q <- diag(length(v_hat_supp[[i]]))
      }else{
        V <- qr.Q(qr(V))
        Q <- diag(length(v_hat_supp[[i]])) - tcrossprod(V)
      }
      if(cov == TRUE){
        QPSPQ <- Q%*%data[v_hat_supp[[i]],v_hat_supp[[i]]]%*%Q
      }else{
        XPQ <- data[,v_hat_supp[[i]]]%*%Q
        QPSPQ <- 1/nrow(XPQ)*crossprod(XPQ)
      }
      eig <- eigen(QPSPQ, TRUE)
      v_hat[v_hat_supp[[i]],i] <- eig$vectors[,1]
      lambda_hat[i] <- eig$values[1]
    }
  }
  
  output <- list("vector" = v_hat, "value" = lambda_hat)
  
  return(output)
}