\name{syn.lognorm, syn.sqrtnorm, syn.cubertnorm}
\alias{syn.lognorm}
\alias{syn.sqrtnorm}
\alias{syn.cubertnorm}
\title{Synthesis by linear regression after transformation of 
       a dependent variable}
\description{
  Generates univariate synthetic data using linear regression
  of an outcome variable transformed by natural logarithm (\code{lognorm}), 
  square root (\code{sqrtnorm}) or cube root (\code{cubertnorm}).
}
\usage{
syn.lognorm(y, x, xp, proper = FALSE, ...) 
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
  fitted linear regression line of transformed \code{y} given \code{x}.
  For proper synthesis first the regression coefficients are drawn
  from normal distribution with mean and variance from the fitted model.
  The synthetic values are transformed back to the original scale.}
\value{A vector of length \code{k} with synthetic values of \code{y}.}
\seealso{
  \code{\link{syn}}, \code{\link{syn.norm}}, \code{\link{syn.normrank}}}
\keyword{datagen}


