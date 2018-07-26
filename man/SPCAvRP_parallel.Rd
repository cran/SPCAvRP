\name{SPCAvRP_parallel}
\alias{SPCAvRP_parallel}
\title{Parallel implementation of the SPCAvRP algorithm}
\description{Computes \code{l}-sparse leading eigenvector of the sample covariance matrix, by parallel selection of \code{A} projections, where each projection is selected from a group of \code{B} random projections of dimension \code{d}.}
\usage{SPCAvRP_parallel(data, cov = FALSE, l, d = 10, A = 300, B = 100,
center_data = TRUE, cluster_type = "PSOCK", cores = 1, machine_names = NULL)}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix.}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{l}{Desired sparsity level in the final estimator (see Details).}
  \item{d}{The dimension of the random projections.}
  \item{A}{Number of projections over which to aggregate.}
  \item{B}{Number of projections in a group from which to select.}
  \item{center_data}{\code{TRUE} if the data matrix should be centered.}
  \item{cluster_type}{Can be \code{"PSOCK"} or \code{"FORK"} (cf. package \code{"parallel"}).}
  \item{cores}{Number of cores to use if \code{clustertype=="FORK"}.}
  \item{machine_names}{Names of computers on the network if \code{clustertype=="PSOCK"}.}
}
\details{This function implements the parallelised SPCAvRP algorithm, by calling \code{'select_projection'} \code{A} times in parallel. We recommend to use this function if \code{p}, \code{A} and \code{B} are large; otherwise use \code{\link{SPCAvRP}}. 

If the true sparsity level \code{k} is known, use \code{d = k} and \code{l = k}.  If \code{k} is unknown, the default choice for \code{d} is \code{10}, while \code{l} can take an array of different values and then the estimators of the corresponding sparsity levels are computed. 

It is desirable to choose \code{A} as big as possible subject to the computational budget. In general, we suggest using \code{A = 300} and \code{B = 100} when the dimension of data is a few hundreds, while \code{A = 600} and \code{B = 200} when the dimension is on order of \code{1000}. 

If \code{center_data == TRUE} and \code{data} is given as a data matrix, the first step is to center it by executing \code{scale(data, center_data, FALSE)}, which subtracts the column means of \code{data} from their corresponding columns.}
\value{Returns a list of two elements:
\item{vector}{A vector or a matrix with \code{length(l)} columns as the estimated eigenvectors of sparsity level \code{l}.}
\item{value}{An array with \code{length(l)} estimated eigenvalues.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\seealso{\code{\link{SPCAvRP}}}
\examples{
p <-  100 # dimension of data
k <- 10 # true sparsity level
n <- 1000 # number of observations
v1 <- c(rep(1/sqrt(k), k), rep(0,p-k)) # leading eigenvector
Sigma <- 2*tcrossprod(v1) + diag(p) # population covariance
mu <- rep(0, p) # population mean
X <- mvrnorm(n, mu, Sigma) # data matrix

spca <- SPCAvRP_parallel(data = X, cov = FALSE, l = k, d = k, A = 200, B = 70, 
                          center_data = FALSE, cluster_type = "PSOCK")
}

