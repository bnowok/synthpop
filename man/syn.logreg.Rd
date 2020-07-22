\name{syn.logreg}
\alias{syn.logreg}
\title{Synthesis by logistic regression}
\description{
  Generates univariate synthetic data for binary or binomial 
  response variable using logistic regression model.
}
\usage{
syn.logreg(y, x, xp, denom = NULL, denomp = NULL, proper = FALSE, ...) 
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{denom}{an original denominator vector of length \code{n} for a binomial 
    regression model.}
  \item{denomp}{a synthesised denominator vector of length \code{k} for 
    a binomial regression model.}
  \item{proper}{a logical value specifying whether proper synthesis 
    should be conducted. See details.}
  \item{\dots}{additional parameters.}
}
\details{
  Synthesis for binary response variables by the non-Bayesian or 
  approximate Bayesian logistic regression model. The non-Bayesian 
  method consists of the following steps: \enumerate{\item Fit a
  logistic regression to the original data. 
  \item Calculate predicted inverse logits for synthesied covariates.
  \item Compare the inverse logits to a random (0,1) deviate and get synthetic 
    values.} 
  The Bayesian version (for proper synthesis) includes additional step before 
  computing inverse logits, namely drawing coefficients from normal distribution 
  with mean and variance estimated in step 1.
  
  The method relies on the standard \code{glm.fit} function.
  Warnings from \code{glm.fit} are suppressed. Perfect prediction
  is handled by the data augmentation method.
}
\value{
  A list with two components:
    \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
    \item{fit}{a summary of the model fitted to the observed data and used
      to produce synthetic values.}
}
\seealso{
  \code{\link{syn}}, \code{\link[stats]{glm}}
}
\keyword{datagen}
