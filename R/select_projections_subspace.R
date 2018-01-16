select_projections_subspace <- function( data            # a list of projected covariances to select from
                                       , rand_ind        # corresponding projections
                                       , s               # the dimension of desired subspace (<=ncol(rand_ind)) 
                                       , p               # original dimension
                                       , A = NULL        # number of projections to be aggregated 
)
  # Output : v_hat_stars # a matrix of dimension (sxp)xA with its columns as c(v_1,...,v_s), where v_r is the eigenvector of 
  #          P_a Sigma_hat P_a and P_a is the projection yielding the largest r-th eigenvalue among B different random projections 
{ 
  N <- nrow(rand_ind)
  if( is.null(A) ){ A <- floor(3*sqrt(N/3)) }
  B <- floor(N/A)
  d <- ncol(rand_ind)
  
  # Selection of good projections
  v_hat_stars <- sapply(1:A, function(a){
    group_decomposition <- sapply(1:B, function(b){
      eig <- eigen(data[[B*(a-1)+b]], TRUE)
      return(list(eig$values[1:s],eig$vectors[,1:s]))})
    lambdas_hat <- matrix(unlist(group_decomposition[1,]),nrow = B, byrow = TRUE)
    b_stars <- max.col(t(lambdas_hat))
    v_hat_star <- rep(0,s*p)
    for(i in 1:s){
      v_hat_star[(i-1)*p+rand_ind[B*(a-1)+b_stars[i],]] <- group_decomposition[2,][[b_stars[i]]][,i]
    }
    
    return(v_hat_star)
  }) 
  
  
  return(v_hat_stars) 
}

