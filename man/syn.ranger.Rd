\name{syn.ranger}
\alias{syn.ranger}
\title{Synthesis with a fast implementation of random forests}
\description{
Generates univariate synthetic data using a fast implementation of 
random forests. It uses \code{\link{ranger}} function 
from the \pkg{ranger} package.
}
\usage{
syn.ranger(y, x, xp, smoothing = "", proper = FALSE, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{smoothing}{smoothing method for continuous variables.}
  \item{proper}{...}
  \item{\dots}{additional parameters passed to 
    \code{\link[ranger]{ranger}}.}    % !!! check and update 
}
\details{...}
\value{
A list with two components:
  \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
  \item{fit}{the model fitted to the observed data that was used to produce 
    synthetic values.}
}
\references{...}
\seealso{
  \code{\link{syn}}, \code{\link{syn.rf}},
  \code{\link{syn.bag}}, \code{\link{syn.cart}},
  \code{\link[ranger]{ranger}}}
\keyword{datagen}

