\name{glasso.miss}
\alias{glasso.miss}
\title{Sparse Covariance Selection by glasso with EM}
\description{A function to achieve covariance selection by glasso with EM, missing data is allowed}
   
\usage{
glasso.miss(Y, emIter, rho, ...)
}

\arguments{
    \item{Y}{numeric matrix. Columns are for variables and rows are for samples.  }
    \item{emIter}{numeric value. The maximum number of iteration for EM} 
	\item{rho}{The l1 norm tuning parameter.}
    \item{...}{additional control parameter passed to \link[=glasso]{glasso}}
}        



\value{
a list of following comopnents.
    \item{Y.imputed}{The sample matrix Y with missing data imputed by EM}
    \item{bic}{the BIC(Bayesian Information Criterion) of current fit.} 
	\item{...}{other values returned by \link[=glasso]{glasso}}
}

\references{

He, S.Y., Wang, X., and Yuan, W.(2012), Discovering Co-movement Structure of Chinese Stock Market by SPACE method with EM

Jerome Friedman, Trevor Hastie and Robert Tibshirani (2007). Sparse inverse covariance estimation with the lasso. Biostatistics 2007. 

Meinshausen, N. and Buhlmann, P.(2006) High dimensional graphs and variable selection with the lasso. Annals of Statistics,34, p1436-1462.

Daniela Witten, Jerome Friedman, and Noah Simon (2011). New insights and faster computations for the graphical lasso. To appear in Journal of Computational and Graphical Statistics.

}

\author{Shiyuan He}

\keyword{methods}

\examples{
data(finStocksCn)  ##data of finance sector of Chinese Stock Market, from January 4th, 2011 to July 26th, 2011
finStocksCn$names ##stock names
y.m<-scale(finStocksCn$returns)
n=nrow(y.m)
p=ncol(y.m)
res=glasso.miss(y.m,rho=0.030,emIter=25,penalize.diagonal=FALSE)
res$bic  ##bic returned
}
