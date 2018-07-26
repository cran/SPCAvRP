SPCAvRP_subspace <- function( data               # either the data matrix or the sample covariance matrix
                            , cov = FALSE        # TRUE if data is given as a sample covariance matrix 
                            , s                  # the dimension of the eigenspace (the number of principal components)
                            , l                  # sparsity level of the eigenvectors (array of lenght s) 
                            , d = 10             # dimension of the random projections, must be >= s (if k is known, use d = k)
                            , A = 300            # number of projections to aggregate 
                            , B = 100            # number of projections in a group from which to select
                            , center_data = TRUE # TRUE if the data matrix should be centered 
)
  # output : a list of two elements
  #             output$vector : a matrix of dimension pxs with its columns as estimated eigenvectors of sparsity level l 
  #             output$value  : an array of dimension s with estimated eigenvalues
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
  
  # Generate random projections
  rand_ind <- matrix(replicate(A*B,sample.int(p,d)), nrow = A*B, byrow = TRUE)
  
  # Project sample covariance
  cov_projections <- project_covariance(data, cov, rand_ind) 
  
  # Selection of good projections
  v_hat_stars <- select_projections_subspace(cov_projections, rand_ind, s, p, A)
  
  # Aggregation of good projections
  w <- NULL; initial_ranking <- NULL; top_coord <- NULL
  for(i in 1:s){
    w[[i]] <- (1/A)*rowSums(abs(v_hat_stars[((i-1)*p+1):(i*p),]))
    initial_ranking[[i]] <- sort(w[[i]], TRUE, index.return = TRUE)$ix 
    top_coord[[i]] <- initial_ranking[[i]][1:min(sum(l[1:i]),p)] # initially keep more coordinates
  }
  
  # Computing orthogonal estimators
  lambda_hat <- rep(0,s)
  v_hat <- matrix(rep(0,s*p), nrow = p)
  if( cov == TRUE ){
    eig <- eigen(data[top_coord[[1]],top_coord[[1]]], TRUE) 
    lambda_hat[1] <- eig$values[1]  
    v_hat[top_coord[[1]],1] <- eig$vectors[,1]
    HSH <- data
    for(i in 2:s){
      H <- diag(length(top_coord[[i-1]]))-tcrossprod(v_hat[top_coord[[i-1]],i-1])
      HSH[top_coord[[i-1]],top_coord[[i-1]]] <- H%*%HSH[top_coord[[i-1]],top_coord[[i-1]]]%*%H
      eigenv_init <- eigen(HSH[top_coord[[i]],top_coord[[i]]], TRUE)$vectors[,1]
      top_coord[[i]] <- top_coord[[i]][which(unlist(sapply(eigenv_init, function(x){ if(!is.element(abs(x), head(sort(abs(eigenv_init), TRUE), l[i]))) {x=0} else {x=1} })) == 1)] # update best coordinates
      V <- NULL
      for(j in seq(i-1,1,-1)){
        u <- v_hat[top_coord[[i]],j]
        if(norm(u,"2")!=0){V <- cbind(V,u/norm(u,"2"))}
      }
      if(is.null(V)){
        Q <- diag(length(top_coord[[i]]))
      }else{
        V <- qr.Q(qr(V))
        Q <- diag(length(top_coord[[i]])) - tcrossprod(V)
      }
      QPSPQ <- Q%*%data[top_coord[[i]],top_coord[[i]]]%*%Q
      eig <- eigen(QPSPQ, TRUE)
      v_hat[top_coord[[i]],i] <- eig$vectors[,1]
      lambda_hat[i] <- eig$values[1]
    }
  }else{
    eig <- eigen(1/nrow(data)*crossprod(data[,top_coord[[1]]]), TRUE)
    lambda_hat[1] <- eig$values[1]  
    v_hat[top_coord[[1]],1] <- eig$vectors[,1]
    XH <- data
    for(i in 2:s){
      H <- diag(length(top_coord[[i-1]]))-tcrossprod(v_hat[top_coord[[i-1]],i-1])
      XH[,top_coord[[i-1]]] <- XH[,top_coord[[i-1]]]%*%H
      eigenv_init <- eigen(1/nrow(XH)*crossprod(XH[,top_coord[[i]]]), TRUE)$vectors[,1]
      top_coord[[i]] <- top_coord[[i]][which(unlist(sapply(eigenv_init, function(x){ if(!is.element(abs(x), head(sort(abs(eigenv_init), TRUE), l[i]))) {x=0} else {x=1} })) == 1)] # update best coordinates
      V <- NULL
      for(j in seq(i-1,1,-1)){
        u <- v_hat[top_coord[[i]],j]
        if(norm(u,"2")!=0){V <- cbind(V,u/norm(u,"2"))}
      }
      if(is.null(V)){
        Q <- diag(length(top_coord[[i]]))
      }else{
        V <- qr.Q(qr(V))
        Q <- diag(length(top_coord[[i]])) - tcrossprod(V)
      }
      XPQ <- data[,top_coord[[i]]]%*%Q
      eig <- eigen(1/nrow(XPQ)*crossprod(XPQ), TRUE)
      v_hat[top_coord[[i]],i] <- eig$vectors[,1]
      lambda_hat[i] <- eig$values[1]
    }
  }
  
  output <- list("vector" = v_hat, "value" = lambda_hat)
  
  return(output)
}