\name{predict.ccn}
\alias{predict.ccn}
\alias{predict.hsvmclusterpathc}
\alias{predict.mcppathc}
\alias{predict.hsvmclusterscadpathc}
\title{make predictions from a "ccn" object.}

\description{
Similar to other predict methods, this functions predicts fitted values and class labels from a fitted \code{\link{gcdclust}} object.
}

\usage{
\method{predict}{ccn}(object, newx, s = NULL,
type=c("class","link"), ...)
}

\arguments{
		\item{object}{fitted \code{\link{ccn}} model object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. NOTE: \code{newx} must be a matrix, \code{predict} function does not accept a vector or other formats of \code{newx}.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model. }
		\item{type}{
		type of prediction required. 
		\itemize{
				\item Type \code{"link"} gives the
				linear predictors for classification problems and gives predicted response for regression problems.
				\item Type \code{"class"} produces the class label corresponding to the maximum probability. Only available for classification problems.}
		}
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


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
