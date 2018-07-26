\name{SPCAvRP_deflation}
\alias{SPCAvRP_deflation}
\title{Computes the leading eigenvectors using the modified deflation scheme}
\description{Computes \code{s} leading eigenvectors of the sample covariance matrix which are sparse and orthogonal, using the modified deflation scheme in conjunction with the SPCAvRP algorithm.}
\usage{SPCAvRP_deflation(data, cov = FALSE, s, l, d = 10, 
A = 300, B = 100, center_data = TRUE)}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix.}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{s}{The number of eigenvectors to compute.}
  \item{l}{The array of lenght \code{s} with the desired sparsity levels in the final estimators.}
  \item{d}{The dimension of the random projections.}
  \item{A}{Number of projections over which to aggregate.}
  \item{B}{Number of projections in a group from which to select.}
  \item{center_data}{\code{TRUE} if the data matrix should be centered.}
}
\details{This function implements the modified deflation scheme in conjunction with the SPCAvRP in order to compute \code{s} sparse eigenvectors that are orthogonal. If possible, use \code{\link{SPCAvRP_subspace}} instead. 
  
If the true sparsity level is known and for each component is equal to \code{k}, use \code{d = k} and \code{l = rep(k,s)}. Sparsity levels of different components may take different values. If \code{k} is unknown, appropriate \code{k} could be chosen from an array of different values by inspecting the explained variance for one component at the time and by using \code{SPCAvRP} in a combination with the deflation scheme implemented in \code{SPCAvRP_deflation}. 
  
It is desirable to choose \code{A} as big as possible subject to the computational budget. In general, we suggest using \code{A = 300} and \code{B = 100} when the dimension of data is a few hundreds, while \code{A = 600} and \code{B = 200} when the dimension is on order of \code{1000}. 
  
If \code{center_data == TRUE} and \code{data} is given as a data matrix, the first step is to center it by executing \code{scale(data, center_data, FALSE)}, which subtracts the column means of \code{data} from their corresponding columns.}
\value{Returns a list of two elements:
  \item{vector}{A matrix whose \code{s} columns are the estimated eigenvectors.}
  \item{value}{An array with \code{s} estimated eigenvalues.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\seealso{\code{\link{SPCAvRP}}, \code{\link{SPCAvRP_subspace}}}
\examples{
p <- 50
k <- 8
theta <- 40
v1 <- c(rep(1/sqrt(k), k), rep(0, p-k))
theta2 <- 20
v2 <- c(rep(0,4), 1/sqrt(k), -1/sqrt(k), 1/sqrt(k), -1/sqrt(k), rep(1/sqrt(k),4), rep(0,p-12))
theta3 <- 5
v3 <- c(rep(0,6), 1/sqrt(k), -rep(1/sqrt(k),4), rep(1/sqrt(k),3), rep(0,p-14))
Sigma <- diag(p) + theta*tcrossprod(v1) + (theta2)*tcrossprod(v2) + (theta3)*tcrossprod(v3) 
mu <- rep(0, p)
n <- 2000
X <- mvrnorm(n, mu, Sigma)

spcarp <- SPCAvRP_deflation(data = X, cov = FALSE, s = 1, l = k, d = k, 
                            A = 300, B = 100, center_data = FALSE)
}

