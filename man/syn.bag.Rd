\name{syn.bag}
\alias{syn.bag}
\title{Synthesis with bagging}
\description{
Generates univariate synthetic data using bagging. It uses 
\code{\link{randomForest}} function from the \pkg{randomForest} package with 
number of sampled predictors equal to number of all predictors. 
}
\usage{
syn.bag(y, x, xp, smoothing = "", proper = FALSE, ntree = 10, ...)
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
  \code{\link{syn}}, \code{\link{syn.rf}}, \code{\link{syn.cart}},
  \code{\link[randomForest]{randomForest}}, \code{\link{syn.smooth}}
}
\keyword{datagen}

