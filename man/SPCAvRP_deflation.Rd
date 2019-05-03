\name{SPCAvRP_deflation}
\alias{SPCAvRP_deflation}
\title{Computes multiple principal components using our modified deflation scheme}
\description{Computes \code{m} leading eigenvectors of the sample covariance matrix which are sparse and orthogonal, using the modified deflation scheme in conjunction with the SPCAvRP algorithm.}
\usage{SPCAvRP_deflation(data, cov = FALSE, m, l, d = 20, 
A = 600, B = 200, center_data = TRUE)}
\arguments{
  \item{data}{Either the data matrix (\code{p x n}) or the sample covariance matrix (\code{p x p}).}
  \item{cov}{\code{TRUE} if data is given as a sample covariance matrix.}
  \item{m}{The number of principal components to estimate.}
  \item{l}{The array of length \code{m} with the desired sparsity of \code{m} principle components (see Details).}
  \item{d}{The dimension of the random projections (see Details).}
  \item{A}{Number of projections over which to aggregate (see Details).}
  \item{B}{Number of projections in a group from which to select (see Details).}
  \item{center_data}{\code{TRUE} if the data matrix should be centered (see Details).}
}
\details{This function implements the modified deflation scheme in conjunction with SPCAvRP (Algorithm 2 in the reference given below). 
  
If the true sparsity level is known and for each component is equal to \code{k}, use \code{d = k} and \code{l = rep(k,m)}. Sparsity levels of different components may take different values. If \code{k} is unknown, appropriate \code{k} could be chosen from an array of different values by inspecting the explained variance for one component at the time and by using \code{SPCAvRP} in a combination with the deflation scheme implemented in \code{SPCAvRP_deflation}. 
  
It is desirable to choose \code{A} (and \code{B = ceiling(A/3)}) as big as possible subject to the computational budget. In general, we suggest using \code{A = 300} and \code{B = 100} when the dimension of data is a few hundreds, while \code{A = 600} and \code{B = 200} when the dimension is on order of \code{1000}. 
  
If \code{center_data == TRUE} and \code{data} is given as a data matrix, the first step is to center it by executing \code{scale(data, center_data, FALSE)}, which subtracts the column means of \code{data} from their corresponding columns.}
\value{Returns a list of two elements:
  \item{vector}{A matrix whose \code{m} columns are the estimated eigenvectors.}
  \item{value}{An array with \code{m} estimated eigenvalues.}
}
\references{Milana Gataric, Tengyao Wang and Richard J. Samworth (2018) Sparse principal component analysis via random projections
\url{https://arxiv.org/abs/1712.05630}}
\author{Milana Gataric, Tengyao Wang and Richard J. Samworth}
\seealso{\code{\link{SPCAvRP}}, \code{\link{SPCAvRP_subspace}}}
\examples{
p <- 50 # data dimension
k <- 8  # true sparsity of each component
v1 <- 1/sqrt(k)*c(rep(1, k), rep(0, p-k)) # first principal compnent (PC)
v2 <- 1/sqrt(k)*c(rep(0,4), 1, -1, 1, -1, rep(1,4), rep(0,p-12)) # 2nd PC
v3 <- 1/sqrt(k)*c(rep(0,6), 1, -rep(1,4), rep(1,3), rep(0,p-14)) # 3rd PC
Sigma <- diag(p) + 40*tcrossprod(v1) + 20*tcrossprod(v2) + 5*tcrossprod(v3) # population covariance 
mu <- rep(0, p) # population mean
n <- 2000 # number of observations
loss = function(u,v){
  sqrt(abs(1-sum(v*u)^2))
}
loss_sub = function(U,V){
  U<-qr.Q(qr(U)); V<-qr.Q(qr(V))
  norm(tcrossprod(U)-tcrossprod(V),"2")
}
set.seed(1)
X <- mvrnorm(n, mu, Sigma) # data matrix

spcavrp.def <- SPCAvRP_deflation(data = X, cov = FALSE, m = 2, l = rep(k,2), 
                                 d = k, A = 200, B = 70, center_data = FALSE)
subspace_estimation<-data.frame(
  loss_sub(matrix(c(v1,v2),ncol=2),spcavrp.def$vector),
  loss(spcavrp.def$vector[,1],v1),
  loss(spcavrp.def$vector[,2],v2),
  crossprod(spcavrp.def$vector[,1],spcavrp.def$vector[,2]))
colnames(subspace_estimation)<-c("loss_sub","loss_v1","loss_v2","inner_prod")
rownames(subspace_estimation)<-c("")
print(subspace_estimation)
}

