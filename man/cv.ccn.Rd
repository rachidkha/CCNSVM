\name{cv.ccn}
\alias{cv.ccn}
\alias{cv.hsvmclusterpathc}
\alias{cv.mcppathc}
\alias{cv.hsvmclusterpathc}
\title{Cross-validation for CCNSVM model}
\description{Does k-fold cross-validation for ccn, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{CCNSVM} package.}
\usage{
cv.ccn(x, y, KK,lambda, pred.loss, nfolds, foldid, delta,...)
}
\arguments{
		\item{x}{\code{x} matrix as in \code{\link{ccn}}.}
		\item{y}{response variable or class label \code{y} as in \code{\link{ccn}}.}
		\item{KK}{the number of clusters determinited by the user}
		\item{lambda}{optional user-supplied lambda sequence; default is
		\code{NULL}, and \code{\link{ccn}} chooses its own sequence.}
		\item{nfolds}{number of folds - default is 5. Although \code{nfolds}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
		\item{foldid}{an optional vector of values between 1 and \code{nfold}
		identifying what fold each observation is in. If supplied,
		\code{nfold} can be missing.}
		\item{pred.loss}{loss function to use for cross-validation error. Valid options are:
		\itemize{
		\item \code{"loss"} Margin based loss function. When use least square loss \code{"ls"}, it gives mean square error (MSE). When use expectile regression loss \code{"er"}, it gives asymmetric mean square error (AMSE).
		\item \code{"misclass"} only gives misclassification error. }
		Default is \code{"loss"}.}
		\item{delta}{parameter \eqn{\delta}{delta} only used in HHSVM loss function for computing margin based loss function, only available for \code{pred.loss = "loss"}.}
		\item{\dots}{other arguments that can be passed to ccn.}

}
\details{The function runs \code{\link{ccn}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the
folds are computed. 
}
\value{an object of class \code{\link{cv.ccn}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting
		purposes).}
		\item{ccn.fit}{a fitted \code{\link{ccn}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
}

\author{Kharoubi, R, Oualkacha, K and Mkhadri, A\cr
Maintainer: Kharoubi Rachid  <kharoubi.rachid@courrier.uqam.ca>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}


\keyword{models}
\keyword{Classification}

