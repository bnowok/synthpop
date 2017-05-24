\name{syn.normrank}
\alias{syn.normrank}
\title{Synthesis by normal linear regression preserving
  the marginal distribution}
\description{
  Generates univariate synthetic data using linear regression analysis
  and preserves the marginal distribution. Regression is carried out on
  Normal deviates of ranks in the original variable. Synthetic values are 
  assigned from the original values based on the synthesised ranks 
  that are transformed from their synthesised Normal deviates.
}
\usage{
syn.normrank(y, x, xp, smoothing, proper = FALSE, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{smoothing}{smoothing method. See details.}
  \item{proper}{a logical value specifying whether proper synthesis 
    should be conducted. See details.}
  \item{\dots}{additional parameters.}
}
\details{
  First generates synthetic values of Normal deviates of ranks of  
  the values in \code{y} using the spread around the fitted 
  linear regression line of Normal deviates of ranks given \code{x}. 
  Then synthetic Normal deviates of ranks are transformed back to 
  get synthetic ranks which are used to assign values from
  \code{y}. 
  For proper synthesis first the regression coefficients 
  are drawn from normal distribution with mean and variance
  from the fitted model.
  A Guassian kernel smoothing can be applied by setting smoothing parameter 
  to \code{"density"}. It is recommended as a tool to decrease the disclosure 
  risk.}
\value{A vector of length \code{k} with synthetic values of \code{y}.}
\seealso{
  \code{\link{syn}}, \code{\link{syn.norm}}, \code{\link{syn.lognorm}}}
\keyword{datagen}
