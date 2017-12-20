select_projection <- function(  data       # either data matrix or sample covariance matrix
                              , cov = TRUE # FALSE if data is given as a data matrix
                              , p          # dimension of the samples
                              , d          # dimension of the random projections
                              , B          # number of random projections to generate and select from
) 
  # Output is v_hat_star : eigenvector of P Sigma_hat P, where P is the projection yielding the largest eigenvalue among B different random projections 
{
  rand_ind <- as.vector(replicate(B,sample.int(p,d))) 
  
  group_decomposition <- sapply(1:B, function(b){
    if(cov == TRUE){
      eig <- eigen(data[rand_ind[d*(b-1)+1:d],rand_ind[d*(b-1)+1:d]], TRUE)
    }else{
      eig <- eigen(1/nrow(data)*crossprod(data[,rand_ind[d*(b-1)+1:d]]), TRUE)
    }
    value <- eig$values[1]
    vector <- eig$vectors[,1]                                 
    return(list(value,vector))})
  lambda_hat <- as.matrix(group_decomposition[1,])
  v_hat <- t(matrix(unlist(group_decomposition[2,]), nrow = B, ncol = d, byrow = TRUE))
  
  b_star <- which.max(lambda_hat)
  v_hat_star <- rep(0,p)
  v_hat_star[rand_ind[d*(b_star-1)+1:d]] <- v_hat[,b_star]
  
  return(v_hat_star)
}
