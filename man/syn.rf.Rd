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
  \item{smoothing}{smoothing method for numeric variable. See 
    \code{\link{syn.smooth}}.}
  \item{proper}{for proper synthesis (\code{proper = TRUE}) a model is fitted 
    to a bootstrapped sample of the original data.}
  \item{ntree}{number of trees to grow.}
  \item{\dots}{additional parameters passed to 
    \code{\link[randomForest]{randomForest}}.}  
}
\details{...}
\value{
A list with two components:
  \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
  \item{fit}{the fitted model which is an object of class \code{randomForest.}}
}
\references{...}
\seealso{
  \code{\link{syn}}, \code{\link{syn.rf}},
  \code{\link{syn.bag}}, \code{\link{syn.cart}},
  \code{\link[randomForest]{randomForest}},
  \code{\link{syn.smooth}}
}
\keyword{datagen}

