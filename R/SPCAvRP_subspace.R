SPCAvRP_subspace <- function( data               # either the data matrix or the sample covariance matrix
                            , cov = FALSE        # TRUE if data is given as a sample covariance matrix
                            , m                  # the dimension of the eigenspace (i.e. the number of principal components)
                            , l                  # sparsity level of the eigenspace (i.e. the number of non-zero rows in output$vector)
                            , d = 20             # dimension of the random projections, must be >= s (if k is known, use d = k)
                            , A = 600            # number of projections to aggregate
                            , B = 200            # number of projections in a group from which to select
                            , center_data = TRUE # TRUE if the data matrix should be centered
)
  # output : a list of three elements
  #             output$vector : a matrix of dimension ncol(data)xm with its columns as the estimated eigenvectors, which has l non-zero rows
  #             output$value  : an array of dimension m with estimated eigenvalues
  #             output$importance_scores  : an array of length ncol(data) with importance scores for each row 1,...,ncol(data)
{
  #########################################
  ################## the required functions:
  
  ## Project the sample covariance matrix along given axis-aligned projections:
  project_covariance <- function(     data          # either data matrix or sample covariance matrix
                                    , cov           # TRUE if data is given as a sample covariance matrix 
                                    , rand_ind      # projections (indices of non-zero elements)
  ) # output: projections - a list of nrow(rand_ind) projected sample covariance matrices of dimension ncol(rand_ind) x ncol(rand_ind)
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
  
  ## Select A projections yielding the largest r-th eigenvalue among B different random projections, for r in c(1:s) :
  select_projections_subspace <- function(   cov_projections # a list of projected covariances to select from
                                           , rand_ind        # corresponding projections
                                           , s               # the dimension of desired subspace (<=ncol(rand_ind)) 
                                           , p               # original dimension
                                           , A = NULL        # number of projections to be aggregated 
  ) # output is v_hat_stars - a matrix of dimension (sxp)xA with its columns as c(v_1,...,v_s), where v_r is the eigenvector of P_a Sigma_hat P_a and P_a is the projection yielding the largest r-th eigenvalue among B different random projections 
  { 
    N <- nrow(rand_ind)
    if( is.null(A) ){ A <- floor(3*sqrt(N/3)) }
    B <- floor(N/A)
    d <- ncol(rand_ind)
    # Selection of good projections
    v_lambda_hat_stars <- sapply(1:A, function(a){
      group_decomposition <- sapply(1:B, function(b){
        eig <- eigen(cov_projections[[B*(a-1)+b]], TRUE)
        if(s == d){eig$values<-append(eig$values,0)}
        return(list(eig$values[1:(s+1)],eig$vectors[,1:s]))})
      lambdas_hat <- matrix(unlist(group_decomposition[1,]),nrow = B, byrow = TRUE)
      b_star <- which.max(colSums(t(lambdas_hat[,1:s])))
      v_hat_star <- rep(0,s*p)
      for(i in 1:s){
        v_hat_star[(i-1)*p+rand_ind[B*(a-1)+b_star,]] <- group_decomposition[2,][[b_star]][,i]
      }
      return(append(v_hat_star,lambdas_hat[b_star,]))
    }) 
    v_hat_stars <- v_lambda_hat_stars[1:(s*p),]
    lambda_hat_stars <- v_lambda_hat_stars[(s*p)+1:(s+1),]
    return(list(v_hat_stars,lambda_hat_stars)) 
  }
  
  ######################################
  ################## the main procedure:
  
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
  SP <- select_projections_subspace(cov_projections, rand_ind, m, p, A)
  v_hat_stars <- SP[[1]]
  lambda_hat_stars <- SP[[2]]

  # Aggregation of good projections by computting the importance scores w
  w<-rep(0,p)
  for(r in 1:m){
      w<-w+rowMeans(matrix(rep(lambda_hat_stars[r,]-lambda_hat_stars[m+1,],p),nrow=p,byrow=TRUE)*v_hat_stars[((r-1)*p+1):(r*p),]^2)
  }
  ranking<-sort(w, TRUE, index.return = TRUE)$ix
  top_coord<-ranking[1:l]

  v_hat <- matrix(rep(0,m*p), nrow = p)
  eig <- eigen(data[top_coord,top_coord], TRUE)
  v_hat[top_coord,] <- eig$vectors[,1:m]
  lambda_hat <- eig$values[1:m]

  output <- list("vector" = v_hat, "value" = lambda_hat, "importance_scores" = w)

  return(output)
}
