\name{syn.pmm}
\alias{syn.pmm}
\title{Synthesis by predictive mean matching}
\description{
  Generates univariate synthetic data using predictive mean
  matching.
}
\usage{
syn.pmm(y, x, xp, proper = FALSE, ...)
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
  Synthesis of \code{y} by predictive mean matching. The procedure
  is as follows: \enumerate{ \item Fit a linear regression 
  to the original data. \item Compute predicted values \code{y.hat}
  and \code{ysyn.hat} for the original \code{x} and synthesised 
  \code{xp} covariates respectively.
  \item For each predicted value \code{ysyn.hat} find donor 
  observations with the closest predicted values \code{y.hat} 
  (ties are broken by random selection), randomly sample one of 
  them and take its observed value \code{y} as the synthetic value. 
  }
  The Bayesian version (for proper synthesis) includes additional 
  step before computing predicted values: \itemize{\item{Draw coefficients 
  from normal distribution with mean and variance
  estimated in step 1 and use them to calculate predicted
  values for the synthesised covariates.}}
}
\value{
A list with two components:
  \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
  \item{fit}{a data frame with regression coefficients and error estimates.}
}
%\references{...}
\seealso{
  \code{\link{syn}}
}
\keyword{datagen}
