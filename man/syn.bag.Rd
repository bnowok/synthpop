\name{syn.bag}
\alias{syn.bag}
\title{Synthesis with bagging}
\description{
Generates univariate synthetic data using bagging. It uses 
\code{\link{randomForest}} function from the \pkg{randomForest} package with 
number of sampled predictors equal to number of all predictors. 
}
\usage{
syn.bag(y, x, xp, smoothing, proper = FALSE, ntree = 10, ...)
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
\value{A vector of length \code{k} with synthetic values of \code{y}.}
\references{...}
\seealso{
  \code{\link{syn}}, \code{\link{syn.rf}}, \code{\link{syn.cart}},
  \code{\link[randomForest]{randomForest}}}
\keyword{datagen}

