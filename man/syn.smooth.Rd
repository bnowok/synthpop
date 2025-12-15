\name{syn.smooth}
\alias{syn.smooth}
\title{syn.smooth}
\description{
  Implements three different smoothing methods for numeric data.
}
\usage{
syn.smooth(ysyn, yobs = NULL, smoothing = "spline", window = 5, ...)
}
\arguments{
  \item{ysyn}{non-missing synthetic data to be smoothed.}
  \item{yobs}{original data used by all methods to determine number of 
    decimal places and by method \code{"density"} to identify 
    top-coding.}
  \item{smoothing}{a character vector that can take values \code{"spline"}, 
    \code{"density"}, or \code{"rmean"}.}
  \item{window}{width of window for running mean.}
  \item{\dots}{additional parameters.}
}
\details{
 Smooths numeric variables by three methods. Default is \code{"spline"} that  
 uses a smoothing spline, others are \code{"density"} that uses a Gaussian 
 kernel density estimator with bandwidth selected using the Sheather-Jones 
 'solve-the-equation' method (see \code{\link[stats]{bw.SJ}}) and \code{"rmean"} 
 that smooths with a running mean of width \code{"window"} (see 
 \code{\link[broman]{runningmean}}). 
}
\value{
 A vector of smoothed values of \code{ysyn}.
}
\seealso{
  \code{\link{syn}}, \code{\link{syn.sample}}, \code{\link{syn.normrank}},
  \code{\link{syn.ctree}}, \code{\link{syn.cart}}, 
  \code{\link{syn.bag}}, \code{\link{syn.rf}}, \code{\link{syn.ranger}},
  \code{\link{syn.nested}}
}
\keyword{smoothing}


