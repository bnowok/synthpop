\name{syn.sample}
\alias{syn.sample}
\title{Synthesis by simple random sampling}
\description{
 Generates a random sample from the observed data.
}
\usage{
syn.sample(y, xp, smoothing, cont.na, proper = FALSE, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{xp}{a target length \code{k} of a synthetic data vector.}
  \item{smoothing}{smoothing method for a continous variable.}
  \item{cont.na}{a vector of codes for missing values for continuous
     variables that should be excluded from smoothing.}  
  \item{proper}{if \code{proper = TRUE} values are sampled from  
    a bootstrapped sample of the original data.}  
  \item{\dots}{additional parameters passed to \code{sample}.}
}
\details{
 A simple random sample with replacement is taken from the 
 observed values in \code{y} and used as synthetic values. 
 A Guassian kernel smoothing can be applied to continuous variables 
 by setting smoothing parameter to \code{"density"}. It is recommended 
 as a tool to decrease the disclosure risk.
}
\value{
 A vector of length \code{k} with synthetic values.
}
\seealso{
  \code{\link{syn}}
}
\keyword{datagen}
