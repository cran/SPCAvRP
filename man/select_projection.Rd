\name{select_projection}
\alias{select_projection}
\title{
Selects the best projection
}
\description{Selects the projection yielding the largest eigenvalue among one group of \code{B} different \code{d}-dimensional axis-aligned projections generated uniformly at random.}
\usage{
select_projection(data, cov = TRUE, p, d, B)
}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix.}
  \item{cov}{\code{FALSE} if data is given as a data matrix.}
  \item{p}{Original dimension of the data.}
  \item{d}{Dimension of the random projections.}
  \item{B}{Number of random projections to generate and select from.}
}
\details{If \code{p} is very large, it might be faster if the data matrix is provided as the input.}
\value{Returns eigenvector \code{v_hat_star} corresponding to the projected covariance yielding the largest eigenvalue among \code{B} different \code{d}-dimensional axis-aligned random projections.}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\examples{
p <-  100
k <- 10
n <- 1000
v1 <- c(rep(1/sqrt(k), k), rep(0,p-k))
Sigma <- 2*tcrossprod(v1) + diag(p)
mu <- rep(0, p)

X <- mvrnorm(n, mu, Sigma)

v_hat_star <- select_projection(data = 1/n*crossprod(X), cov = TRUE, p, d = k, B = 100)
}
