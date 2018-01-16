\name{select_projections_subspace}
\alias{select_projections_subspace}
\title{Selects the best projections for the subspace estimation}
\description{Selects \code{A} projections yielding the largest \code{r}-th eigenvalue among \code{B} different random projections, for \code{r=1:s}.}
\usage{select_projections_subspace(data, rand_ind, s, p, A)}
\arguments{
  \item{data}{A list of projected covariances to select from (in the form of the output of \code{\link{project_covariance}}).}
  \item{rand_ind}{Corresponding projections used to generate \code{data}.}
  \item{s}{The number of eigenvalues to estimate.}
  \item{p}{The original dimension of samples.}
  \item{A}{The number of projections to select.}
}
\value{Returns matrix \code{v_hat_stars} with \code{A} columns that correspond to \code{s} eigenvectors yielding the largest eigenvalues.}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2017) Sparse principal component analysis via random projections
  \url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
