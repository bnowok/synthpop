\name{syn.norm}
\alias{syn.norm}
\title{Synthesis by linear regression}
\description{
  Generates univariate synthetic data using linear regression
  analysis.
}
\usage{
syn.norm(y, x, xp, proper = FALSE, ...) 
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{proper}{a logical value specifying whether proper synthesis 
    should be conducted. See details.}
  \item{\dots}{additional parameters.}
}
\details{
  Generates synthetic values using the spread around the
  fitted linear regression line of \code{y} given \code{x}.
  For proper synthesis first the regression coefficients 
  are drawn from normal distribution with mean and variance
  from the fitted model.
}
\value{
  A list with two components:
   \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
   \item{fit}{a data frame with regression coefficients and error estimates.}
}
\seealso{
  \code{\link{syn}}, \code{\link{syn.normrank}}, \code{\link{syn.lognorm}}}
\keyword{datagen}
