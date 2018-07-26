\name{final_estimator}
\alias{final_estimator}
\title{Computes the leading eigenvector from its support}
\description{Computes the leading eigenvector of the sample covariance matrix given the indices of variables ranked by their importance and desired sparsity level. }
\usage{final_estimator(data, cov, l, ranking)}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix.}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{l}{Desired sparsity of the final estimator (see Details).}
  \item{ranking}{Original variables ranked by their importance.}
}
\details{If true sparsity level \code{k} is known use \code{l = k}. If \code{k} is unknown, \code{l} can be an array of different values and then the eigenvectors of the corrsponding sparsity levels are returned.}
\value{Returns a list of two elements:
\item{vector}{A vector or a matrix with \code{length(l)} columns as the estimated eigenvectors of sparsity level \code{l}.}
\item{value}{An array with \code{length(l)} estimated eigenvalues.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\examples{
p <- 80
k <- 8
n <- 1000
v1 <- c(rep(1/sqrt(k), k), rep(0,p-k))
Sigma <- 2*tcrossprod(v1) + diag(p)
mu <- rep(0, p)
X <- mvrnorm(n, mu, Sigma)
Sigma_hat <- 1/n*crossprod(X)

A <- 200
B <- 100
d <- 10

rand_ind <- matrix(replicate(A*B,sample.int(p,d)), nrow = A*B, byrow = TRUE)
cov_projections <- project_covariance(data = Sigma_hat, cov = TRUE, rand_ind)
ranking <- SPCAvRP_ranking(cov_projections, rand_ind, p, A)

output <- final_estimator(data = Sigma_hat, cov = TRUE, l = (5:11), ranking)
df <- data.frame(5:11,output$value); colnames(df) <- c('l','eigenvalue')
print(df)
}

