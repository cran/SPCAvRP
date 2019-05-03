\name{SPCAvRP}
\alias{SPCAvRP}
\title{Computes the leading eigenvector using the SPCAvRP algorithm}
\description{Computes \code{l}-sparse leading eigenvector of the sample covariance matrix, using \code{A x B} random axis-aligned projections of dimension \code{d}. For the multiple component estimation use \code{\link{SPCAvRP_subspace}} or \code{\link{SPCAvRP_deflation}}.}
\usage{SPCAvRP(data, cov = FALSE, l, d = 20, A = 600, B = 200, 
center_data = TRUE, parallel = FALSE, 
cluster_type = "PSOCK", cores = 1, machine_names = NULL)}
\arguments{
  \item{data}{Either the data matrix (\code{p x n}) or the sample covariance matrix (\code{p x p}).}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{l}{Desired sparsity level in the final estimator (see Details).}
  \item{d}{The dimension of the random projections (see Details).}
  \item{A}{Number of projections over which to aggregate (see Details).}
  \item{B}{Number of projections in a group from which to select (see Details).}
  \item{center_data}{\code{TRUE} if the data matrix should be centered (see Details).}
  \item{parallel}{\code{TRUE} if the selection step should be computed in parallel by uses package \code{"parallel"}.}
  \item{cluster_type}{If \code{parallel == TRUE}, this can be \code{"PSOCK"} or \code{"FORK"} (cf. package \code{"parallel"}).}
  \item{cores}{If \code{parallel == TRUE} and \code{cluster_type == "FORK"}, number of cores to use.}
  \item{machine_names}{If \code{parallel == TRUE}, the names of the computers on the network.}
}
\details{This function implements the SPCAvRP algorithm for the principal component estimation (Algorithm 1 in the reference given below). 

If the true sparsity level \code{k} is known, use \code{l = k} and \code{d = k}. 

If the true sparsity level \code{k} is unknown, \code{l} can take an array of different values and then the estimators of the corresponding sparsity levels are computed. The final choice of \code{l} can then be done by the user via inspecting the explained variance computed in the output \code{value} or via inspecting the output \code{importance_scores}. The default choice for \code{d} is \code{20}, but we suggest choosing \code{d} equal to or slightly larger than \code{l}.  

It is desirable to choose \code{A} (and \code{B = ceiling(A/3)}) as big as possible subject to the computational budget. In general, we suggest using \code{A = 300} and \code{B = 100} when the dimension of data is a few hundreds, while \code{A = 600} and \code{B = 200} when the dimension is on order of \code{1000}. 

If \code{center_data == TRUE} and \code{data} is given as a data matrix, the first step is to center it by executing \code{scale(data, center_data, FALSE)}, which subtracts the column means of \code{data} from their corresponding columns.

If \code{parallel == TRUE}, the parallelised SPCAvRP algorithm is used. We recommend to use this option if \code{p}, \code{A} and \code{B} are very large. 
}
\value{Returns a list of three elements:
\item{vector}{A matrix of dimension \code{p x length(l)} with columns as the estimated eigenvectors of sparsity level \code{l}.}
\item{value}{An array with \code{length(l)} eigenvalues corresponding to the estimated eigenvectors returned in \code{vector}.}
\item{importance_scores}{An array of length p with importance scores for each variable 1 to p.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\examples{
p <- 100  # data dimension
k <- 10   # true sparsity level
n <- 1000 # number of observations
v1 <- c(rep(1/sqrt(k), k), rep(0,p-k)) # true principal component
Sigma <- 2*tcrossprod(v1) + diag(p)    # population covariance
mu <- rep(0, p)                        # population mean
loss = function(u,v){ 
  # the loss function
  sqrt(abs(1-sum(v*u)^2))
}
set.seed(1)
X <- mvrnorm(n, mu, Sigma) # data matrix

spcavrp <- SPCAvRP(data = X, cov = FALSE, l = k, d = k, A = 200, B = 70)
spcavrp.loss <- loss(v1,spcavrp$vector)
print(paste0("estimation loss when l=d=k=10, A=200, B=70: ", spcavrp.loss))

##choosing sparsity level l if k unknown:
#spcavrp.choosel <- SPCAvRP(data = X, cov = FALSE, l = c(1:30), d = 15, A = 200, B = 70)
#plot(1:p,spcavrp.choosel$importance_scores,xlab='variable',ylab='w',
#     main='choosing l when k unknown: \n importance scores w')
#plot(1:30,spcavrp.choosel$value,xlab='l',ylab='Var_l',
#     main='choosing l when k unknown: \n explained variance Var_l')
}

