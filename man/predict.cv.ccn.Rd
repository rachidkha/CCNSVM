\name{predict.cv.ccn}
\alias{predict.cv.ccn}
\title{make predictions from a "cv.ccn" object.}
\description{
This function makes predictions from a cross-validated CCNSVM model,
using the stored \code{"ccn.fit"} object, and the optimal value
chosen for \code{lambda}.
}
\usage{
\method{predict}{cv.ccn}(object, newx, s=c("lambda.1se","lambda.min"),...)
}
\arguments{
		\item{object}{fitted \code{\link{cv.ccn}} object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix. See documentation for \code{predict.ccn}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the value \code{s="lambda.1se"} stored
		on the CV object. Alternatively \code{s="lambda.min"} can be
		used. If \code{s} is numeric, it is taken as the value(s) of
		\code{lambda} to be used.}
		\item{\dots}{not used. Other arguments to predict. } }
\details{This function makes it easier to use the results of
  cross-validation to make a prediction.}
\value{The object returned depends the \dots argument which is passed on
to the \code{\link{predict}} method for \code{\link{ccn}} objects.}

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
