\name{syn.rf}
\alias{syn.rf}
\title{Synthesis with random forest}
\description{
Generates univariate synthetic data using Breiman's random forest algorithm 
classification and regression. It uses \code{\link{randomForest}} function 
from the \pkg{randomForest} package.
}
\usage{
syn.rf(y, x, xp, smoothing = "", proper = FALSE, ntree = 10, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{smoothing}{smoothing method for continuous variables.}
  \item{proper}{...}
  \item{ntree}{number of trees to grow.}
  \item{\dots}{additional parameters passed to 
    \code{\link[randomForest]{randomForest}}.}    % !!! check and update 
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
  \code{\link{syn}}, \code{\link{syn.bag}}, \code{\link{syn.cart}},
  \code{\link[randomForest]{randomForest}}}
\keyword{datagen}

