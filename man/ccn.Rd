\name{ccn}
\alias{ccn}
\title{Fits the regularization paths for large margin classifiers}
\description{Fits a regularization path for large margin classifiers at a sequence of regularization parameters lambda.}
\usage{
ccn(x, y, KK, nlambda = 25, 
		method = c("ccnL1", 
    "ccnscad", "ccnmcp"),
		lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), 
		lambda = NULL, lambda2 = 0, 
		pf = rep(1, nvars), pf2 = rep(1, nvars), exclude, 
		dfmax = nvars + 1, pmax = min(dfmax * 1.2, 
	    nvars), standardize = TRUE, eps = 1e-5, maxit = 700, 
	    delta = 2)
}
\arguments{
		\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}

		\item{y}{response variable. This argument should be a two-level factor for classification.}
		\item{nlambda}{the number of \code{lambda} values - default is 25.}
		\item{method}{a character string specifying the loss function to use, valid options are:
		\itemize{
		\item \code{"ccnL1"} CCN-SVM + L1,
		\item \code{"ccnscad"} CCN-SVM+SCAD,
		\item \code{"ccnmcp"} CCN-SVM + MCP,
		} Default is \code{"hhsvmcluster"}.}
		\item{lambda.factor}{The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p}{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
		close to zero.  If \eqn{N<p}, the default is \code{0.01}.
		A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence.} 
		\item{lambda}{a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
		the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.factor}. Supplying a value of
		\code{lambda} overrides this. It is better to supply
		a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

		\item{lambda2}{regularization parameter \eqn{\lambda_2}{lambda2} for the quadratic penalty of the 
		coefficients.}

		\item{pf}{L1 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. Separate L1 penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential L1 shrinkage. Can be 0 for some variables, which implies
		no L1 shrinkage, and results in that variable always being included in the
		model. Default is 1 for all variables (and implicitly infinity for
		variables listed in \code{exclude}).}
		
		\item{pf2}{L2 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. Separate L2 penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential L2 shrinkage. Can be 0 for some variables, which implies
		no L2 shrinkage. Default is 1 for all variables.}

		\item{exclude}{indices of variables to be excluded from the
		model. Default is none. Equivalent to an infinite penalty factor.}

		\item{dfmax}{limit the maximum number of variables in the
		model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.}

		\item{pmax}{limit the maximum number of variables ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters model through the path, it will be counted only once. Default is \code{min(dfmax*1.2,p)}.}

		\item{standardize}{logical flag for variable standardization, prior to
		fitting the model sequence. If \code{TRUE}, \code{x} matrix is normalized such that \code{x} is centered (i.e. \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}), and sum squares of each column \eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N=1}. If \code{x} matrix is standardized, the ending coefficients will be transformed back to the original scale. Default is \code{FALSE}.}

		\item{eps}{convergence threshold for coordinate majorization descent.  Defaults value is \code{1e-5}.}

		\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda value. Default is 700. If models do not converge, consider increasing \code{maxit}.}

		\item{delta}{the parameter \eqn{\delta}{delta} in the CCNSVM model. The value must be greater than 0. Default is 2.}
		\item{KK}{the number of clusters}
		
}


\value{
An object with S3 class \code{\link{ccn}}.
		\item{call}{the call that produced this object}
		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{dim}{dimension of coefficient matrix (ices)}
		\item{npasses}{total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{Kharoubi, R, Oualkacha, K and Mkhadri, A\cr
Maintainer: Kharoubi Rachid  <kharoubi.rachid@courrier.uqam.ca>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr

}

\keyword{models}
\keyword{Classification}
