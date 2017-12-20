SPCAvRP_parallel <- function(  data                         # either data matrix or sample covariance matrix
                             , cov = FALSE                  # TRUE if data is given as a sample covariance matrix
                             , l                            # sparsity of the final estimator
                                                            # (if k is unknown, this can be an array of at most p different values and then estimators of corrsponding sparsity levels are returned)
                             , d = 10                       # dimension of the random projections (if k is known, use d = k)
                             , A = 300                      # number of groups of projections
                             , B = 100                      # number of projections in a group
                             , datascaling = TRUE           # TRUE if function scale() should be applied to the data matrix with its arguments center and scale 
                             , center = TRUE                # cf. arguments of function scale()
                             , scale = TRUE                 # cf. arguments of function scale()
                             , cluster_type = "PSOCK"       # can be "PSOCK" or "FORK" (cf. package "parallel")
                             , cores =  1
                             , machine_names = NULL         #names of computers on network to form cluster
)
  # output : a list of two elements
  #             output$vector : a vector / matrix with its columns as estimated eigenvectors of different sparsity level l
  #             output$value  : an array of lenght(l) with estimated eigenvalues
{
  if( !cov & datascaling ){
    data <- scale(data, center, scale)
  }
  if(cluster_type == "FORK"){
    cluster = makeCluster(cores, type="FORK")
    #clusterEvalQ(cl = cluster, tools::psnice(value = 19))
  }
  if(cluster_type == "PSOCK"){
    if(is.null(machine_names)){
      cluster <- makePSOCKcluster(1)
    }else{
      cluster <- makePSOCKcluster(names = machine_names)
    }
    clusterExport(cluster, c("select_projection"), envir = environment())
    #clusterEvalQ(cl = cluster, tools::psnice(value = 19))
  }

  p <- ncol(data)
  p_cond <- (p <= 2*d*sqrt(A*B))
  if( p_cond & cov == FALSE ){
    data <- 1/nrow(data)*crossprod(data)
    cov <- TRUE
  }else if( !p_cond & cov == FALSE ){
    print('SPCAvRP: since p > 2*d*sqrt(A*B), computation of p x p sample covariance is avoided')
  }

  # Selection of good projections in parallel
  v_hat_stars <- parSapply(cl = cluster, 1:A, function(a){return(select_projection(data, cov, p, d, B))})
  stopCluster(cluster)

  # Aggregation of good projections
  w <- rowSums(abs(v_hat_stars))/A
  ranking <- sort(w, TRUE, index.return=TRUE)$ix

  # Selection of top l coordinates and compution of the final estimator
  output <- final_estimator(data, cov, l, ranking)

  return(output)
}
