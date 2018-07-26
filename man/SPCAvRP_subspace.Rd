\name{SPCAvRP_subspace}
\alias{SPCAvRP_subspace}
\title{Computes the leading eigenspace using 
the SPCAvRP algorithm for eigenspace estimation}
\description{Computes \code{s} leading eigenvectors of the sample covariance matrix which are sparse and orthogonal, using \code{A} groups of \code{B} random axis-aligned projections of dimension \code{d}.}
\usage{SPCAvRP_subspace(data, cov = FALSE, s, l, d = 10, 
A = 300, B = 100, center_data = TRUE)}
\arguments{
  \item{data}{Either the data matrix or the sample covariance matrix.}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{s}{The dimension of the eigenspace, i.e the number of principal components.}
  \item{l}{The array of length \code{s} with the desired sparsity levels in the final estimators.}
  \item{d}{The dimension of the random projections.}
  \item{A}{Number of projections over which to aggregate.}
  \item{B}{Number of projections in a group from which to select.}
  \item{center_data}{\code{TRUE} if the data matrix should be centered.}
}
\details{This function implements the SPCAvRP algorithm for eigenspace estimation. 
  
If the true sparsity level is known and for each component is equal to \code{k}, use \code{d = k} and \code{l = rep(k,s)}. Sparsity levels of different components may take different values. If \code{k} is unknown, appropriate \code{k} could be chosen from an array of different values by inspecting the explained variance for one component at the time and by using \code{\link{SPCAvRP}} in a combination with \code{\link{SPCAvRP_deflation}}. 
  
It is desirable to choose \code{A} as big as possible subject to the computational budget. In general, we suggest using \code{A = 300} and \code{B = 100} when the dimension of data is a few hundreds, while \code{A = 600} and \code{B = 200} when the dimension is on order of \code{1000}. 

If \code{center_data == TRUE} and \code{data} is given as a data matrix, the first step is to center it by executing \code{scale(data, center_data, FALSE)}, which subtracts the column means of \code{data} from their corresponding columns.}
\value{Returns a list of two elements:
  \item{vector}{A matrix whose \code{s} columns are the estimated eigenvectors.}
  \item{value}{An array with \code{s} estimated eigenvalues.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\seealso{\code{\link{SPCAvRP}}, \code{\link{SPCAvRP_deflation}}}
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

loss = function(u,v){
  sqrt(abs(1-sum(v*u)^2))
}
loss_sub = function(U,V){
  V<-qr.Q(qr(V))
  norm(tcrossprod(U)-tcrossprod(V),"2")
}

s <- 2
spcarp <- SPCAvRP_subspace(data = X, cov = FALSE, s, l = rep(k,s), d = k,
                           A = 200, B = 100, center_data = FALSE)

subspace_estimation<-data.frame(
  loss(spcarp$vector[,1],v1),
  loss(spcarp$vector[,2],v2),
  loss_sub(matrix(c(v1,v2),ncol=s),spcarp$vector),
  crossprod(spcarp$vector[,1],spcarp$vector[,2]))
colnames(subspace_estimation)<-c("loss_1","loss_2","loss_sub","inp")
rownames(subspace_estimation)<-c("")
subspace_estimation
}

