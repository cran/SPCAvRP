final_estimator <- function(  data            # either data matrix or sample covariance matrix
                            , cov             # TRUE if data is given as a sample covariance matrix 
                            , l               # sparsity of the final estimator 
                                              # (if k is unknown, this can be an array of at most p different values and then estimators of corrsponding sparsity levels are returned)
                            , ranking         # p varaibles ranked by their importance
)
  # output : a list of two elements
  #             output$vector : a vector / matrix with its columns as estimated eigenvectors of different sparsity level l 
  #             output$value  : an array of lenght(l) with estimated eigenvalues
{ 
  p <- ncol(data)
  
  l_loop <-sapply(1:length(l), function(i){
    top_coord <- ranking[1:l[i]]
    if( cov == TRUE ){
      eig <- eigen(data[top_coord,top_coord], TRUE)
    }else{
      eig <- eigen(1/nrow(data)*crossprod(data[,top_coord]), TRUE)
    }
    lambda_hat <- eig$values[1]
    v_hat <- rep(0,p)  
    v_hat[top_coord] <- eig$vectors[,1] 
    return(list(v_hat,lambda_hat))
  })
  v_hat <- matrix(unlist(l_loop[1,]), nrow = p, ncol = length(l) , byrow = FALSE)
  if(length(l)==1){v_hat <- as.vector(v_hat)}
  lambda_hat <- unlist(l_loop[2,])
  
  output <- list("vector" = v_hat, "value" = lambda_hat)
  
  return(output) 
}