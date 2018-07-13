\name{print.ccn}
\alias{print.ccn}
\title{print a ccn object}
\description{
Print a summary of the gcdnet path at each step along the path.
}
\usage{
\method{print}{ccn}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{x}{fitted \code{\link{ccn}} object}
  \item{digits}{significant digits in printout}
  \item{\dots}{additional print arguments}
}
\details{
The call that produced the \code{\link{ccn}} object is printed, followed by a two-column matrix with columns \code{Df} and \code{Lambda}. The \code{Df} column is the number of nonzero coefficients.
}
\value{
a two-column matrix, the first columns is the number of nonzero coefficients and the second column is \code{Lambda}.}

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

