\name{SPCAvRP}
\alias{SPCAvRP}
\title{Computes the leading eigenvector using the SPCAvRP algorithm}
\description{Computes \code{l}-sparse leading eigenvector of the sample covariance matrix, using \code{A} groups of \code{B} random axis-aligned projections of dimension \code{d}.}
\usage{SPCAvRP(data, cov = FALSE, l, d = 10, A = 300, B = 100, 
datascaling = TRUE, center = TRUE, scale = TRUE)}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix.}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{l}{Desired sparsity level in the final estimator (see Details).}
  \item{d}{The dimension of the random projections.}
  \item{A}{Number of projections over which to select and aggregate.}
  \item{B}{Number of projections in a group from which to select.}
  \item{datascaling}{\code{TRUE} if function \code{scale()} should be applied to the data matrix with its arguments \code{center} and \code{scale}.}
  \item{center}{If \code{datascaling == TRUE}, this is the argument of function \code{scale()}.}
  \item{scale}{If \code{datascaling == TRUE}, this is the argument of function \code{scale()}.}
}
\details{This function implements the SPCAvRP algorithm. 

If the true sparsity level \code{k} is known, use \code{d = k} and \code{l = k}. If \code{k} is unknown, the default choice for \code{d} is \code{10}, while \code{l} can take an array of different values and then the estimators of the corresponding sparsity levels are computed.

We suggest using \code{A = 300} and \code{B = 100} when the dimension of data is a few hundreds, while \code{A = 600} and \code{B = 200} when the dimension is on order of \code{1000}. 

If \code{datascaling == TRUE} and \code{data} is given as a data matrix, the first step is to scale it by executing \code{scale(data, center, scale)}. By default \code{center == TRUE}, which means that centering is done by subtracting the column means of \code{data} from their corresponding columns; also \code{scale == TRUE}, which means that the scaling is done by dividing the (centered) columns of \code{data} by their standard deviations in case \code{center == TRUE}, and the root mean square otherwise.}
\value{Returns a list of two elements:
\item{vector}{A vector or a matrix with \code{length(l)} columns as the estimated eigenvectors of sparsity level \code{l}.}
\item{value}{An array with \code{length(l)} estimated eigenvalues.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2017) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\seealso{\code{\link{SPCAvRP_parallel}}, \code{\link{SPCAvRP_ranking}}}
\examples{
p <-  100 # dimension of data
k <- 10 # true sparsity level
n <- 1000 # number of observations
v1 <- c(rep(1/sqrt(k), k), rep(0,p-k)) # leading eigenvector
Sigma <- 2*tcrossprod(v1) + diag(p) # population covariance
mu <- rep(0, p) # population mean
X <- mvrnorm(n, mu, Sigma) # data matrix

spca <- SPCAvRP(data = X, cov = FALSE, l = k, d = k, A = 200, B = 70, datascaling = FALSE)
spca$vector
spca$value
}
