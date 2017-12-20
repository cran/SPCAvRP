project_covariance <- function( data          # either data matrix or sample covariance matrix
                              , cov           # TRUE if data is given as a sample covariance matrix 
                              , rand_ind      # projections (indices of non-zero elements)
) 
  # Output : projections # list of nrow(rand_ind) projected sample covariance matrices of dimension ncol(rand_ind) x ncol(rand_ind)
{
  M <- nrow(rand_ind) 
  
  if( cov == TRUE ){
    projections <- sapply(1:M, function(i){
      proj_ind <- rand_ind[i,]
      cov_proj <- data[proj_ind,proj_ind]
      return(list(cov_proj))
    })
  }else{
    projections <- sapply(1:M, function(i){
      proj_ind <- rand_ind[i,]
      cov_proj <- 1/nrow(data)*crossprod(data[,proj_ind])
      return(list(cov_proj))
    })
  }
  
  return(projections)
}