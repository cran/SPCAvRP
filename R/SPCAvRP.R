SPCAvRP <- function( data               # either data matrix or sample covariance matrix
                   , cov = FALSE        # TRUE if data is given as a sample covariance matrix 
                   , l                  # sparsity level of the final estimator 
                                        # (if k is unknown, this can be an array of at most p different values and then estimators of corrsponding sparsity levels are returned)
                   , d = 10             # dimension of the random projections (if k is known, use d = k)
                   , A = 300            # number of groups of projections 
                   , B = 100            # number of projections in a group
                   , center_data = TRUE # TRUE if the data matrix should be centered 
)
  # output : a list of two elements
  #             output$vector : a vector / matrix with its columns as estimated eigenvectors of sparsity level l 
  #             output$value  : an array of lenght(l) with estimated eigenvalues
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
  
  # Ranking of variables
  ranking <- SPCAvRP_ranking(cov_projections, rand_ind, p, A) 
  
  output <- final_estimator(data, cov, l, ranking)
  
  return(output)
}