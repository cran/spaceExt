\name{space.miss}
\alias{space.miss}
\title{Sparse Covariance Selection by SPACE with EM}
\description{A function to estimate partial correlations using SPACE method with EM, missing data is allowed}
   
\usage{
space.miss(Y.m,lam1, lam2=0, sig=NULL, weight=NULL,iter=2, emIter=5,n_iter=1000,t0=0,r=0)
}

\arguments{
\item{Y.m}{numeric matrix. Columns are for variables and rows are for samples. It's recommended to first standardize each column to have mean 0 and \eqn{l_2} norm 1.}
\item{lam1}{numeric value. This is the \eqn{l_1} norm penalty parameter. If the columns of Y.m have norm one,
then the suggested range of lam1 is \eqn{O(n^{3/2}\Phi^{-1}(1-\alpha/(2p^2)))} for small \eqn{\alpha} such as 0.1.}
\item{lam2}{numeric value. If not specified, lasso regression is used in the Joint Sparse Regression Model (JSRM). 
Otherwise, elastic net regression is used in JSRM and \var{lam2} serves as the \eqn{l_2} norm penalty parameter.}
\item{sig}{numeric vector. Its length should be the same as the number of columns of Y.m. It is the vector of \eqn{\sigma^{ii}} (the diagonal of the inverse covariance matrix).
If not specified, \eqn{\sigma^{ii}} will be estimated during the model fitting with initial values rep(1,p).
The number of the iteration of the model fitting (iter) will then be at least 2. Note, the scale of sig does not matter.}
\item{weight}{numeric value or vector. It specifies the weights or the type of weights used for each regression in JSRM. The default value is NULL,
which means all regressions will be weighted equally in the joint model. 
If weight=1, residue variances will be used for weights. 
If weight=2, the estimated degree of each variable will be used for weights.
Otherwise, it should be a positive numeric vector, whose length is equal to the number of columns of Y.m}
\item{iter}{integer. It is the total number of interactions in JSRM for estimating \eqn{\sigma^{ii}} and partial correlations. When sig=NULL and/or weight=NULL or 2, iter should be at least 2.}
\item{emIter}{integer. the maximum number of EM iteration allowed}
\item{n_iter}{integer. the maximum number of interations in JSRM.}
\item{t0}{integer. \eqn{1<=t0<=n}. the time point at which to perform local smoothing estimation.}
\item{r}{positive value. local smoothing parameter. If r=0, then no smoothing is incorporated.}
}

\details{
\code{space.miss} Based on the work of J. Peng, et al(2007), this function allows computing with missing data by 
	implementing EM algorithm. 
	\deqn{
	\frac{1}{2}\sum_{t=1}^n K(r*|t-t_0|)
	[\sum_{i=1}^p w_i(y_{it}-\sum_{j\neq i}\sqrt{\frac{\sigma_{jj}}{\sigma_{ii}}}\rho_{ij}y_{jt})^2 ]
	+\lambda\sum_{i<j} |\rho_ij|
	}
	where  \eqn{K(\cdot)} is the smoothing kernal, \eqn{K(x)=\exp(-x)}.    
}

\value{
	a list of following comopnents.
    \item{Y.imputed}{The sample matrix Y with missing data imputed by EM}
    \item{ParCor}{the estimated partial correlation matrix.}
    \item{sig.fit}{numeric vector of the estimated diagonal \eqn{\sigma^{ii}}.} 
	\item{bic}{BIC for the current estimate.}
}

\references{

He, S.Y., Wang, X., and Yuan, W.(2012), Discovering Co-movement Structure of Chinese Stock Market by SPACE method with EM.

J. Peng, P. Wang, N. Zhou, J. Zhu (2007), Partial Correlation Estimation by Joint Sparse Regression Model. 

Meinshausen, N., and Buhlmann, P. (2006), High Dimensional Graphs and Variable Selection with the Lasso,  Annals of Statistics, 34, 1436-1462.

}

\author{Shiyuan He}

\keyword{methods}

\examples{
data(finStocksCn)  ##data of finance sector of Chinese Stock Market, from January 4th, 2011 to July 26th, 2011
finStocksCn$names ##stock names
y.m<-scale(finStocksCn$returns)
n=nrow(y.m)
p=ncol(y.m)
alpha=0.1
l1=1/sqrt(n)*qnorm(1-alpha/(2*p^2))
res=space.miss(Y.m=y.m,lam1=l1*25,emIter=25,iter=2,weight=2)
res$bic  ##bic returned
}
