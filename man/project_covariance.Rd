\name{project_covariance}
\alias{project_covariance}
\title{Projects the sample covariance}
\description{Projects the sample covariance matrix along given axis-aligned projections.}
\usage{project_covariance(data, cov, rand_ind)}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix (see Details).}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{rand_ind}{Matrix whose rows are the indices of non-zero entries of axis-aligned projections.}
}
\details{If the dimension of data is very large, it might be faster if the data matrix is provided as the input.}
\value{Returns a list of the sample covariance projections:
\item{projections[[1]]}{projected sample covariance along indices of \code{rand_ind[1,]}}
...
\item{projections[[nrow(rand_ind)]]}{projected sample covariance along indices of \code{rand_ind[nrow(rand_ind),]}}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2017) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\examples{
p <-  50 # dimension of data
k <- 5 # sparsity level
n <- 1000 # number of observations
v1 <- c(rep(1/sqrt(k), k), rep(0,p-k))
Sigma_hat <- 1/n*crossprod(mvrnorm(n, rep(0,p), tcrossprod(v1)+diag(p)))

N <- 1000 # number of projections
d <- k # dimension of projections
rand_ind <- matrix(replicate(N,sample.int(p,d)), nrow = N, byrow = TRUE) # axis-aligned projections

cov_projections <- project_covariance(data = Sigma_hat, cov = TRUE, rand_ind)
}
