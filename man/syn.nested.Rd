\name{syn.nested}
\alias{syn.nested}
\title{Synthesis for a variable nested within another variable.}
\description{Synthesizes one variable (\code{y}) from another one (\code{x}) 
  when \code{y} is nested in the categories of \code{x}. A bootstrap 
  sample is created from the original values of \code{y} within each category 
  of \code{xp} (the syntheisied values of the grouping variable).}
\usage{
syn.nested(y, x, xp, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n} for the nested variable.}
  \item{x}{an original data vector of length \code{n} for the variable within 
    which \code{y} is nested.}
  \item{xp}{a vector of length \code{k} with synthetic values of \code{x}.}
  \item{\dots}{additional parameters.}
}
\details{An example would be when \code{x} is a detailed classification 
  of occupations and \code{y} is more detailed sub-classification. It is 
  intended that both \code{x} and \code{y} are categorical (factor) variables. 
  %At present it will not work correctly for numerical variables if there are missing data.
  A warning will be issued if the original y is not nested within x.}
\value{A vector of length \code{k} with synthetic values of \code{y}.}
\keyword{datagen}
