\name{compare}
\alias{compare}
\title{
 Comparison of synthesised and observed data
}
\description{
 A generic function for comparison of synthesised and observed data. 
 The function invokes particular methods which depend on the class 
 of the first argument. 
}
\usage{
 compare(object, data, ...)
}
\arguments{
  \item{object}{a synthetic data object of class \code{synds} 
    or \code{fit.synds}.}
  \item{data}{an original observed data set.}
  \item{\dots}{additional arguments specific to a method.}
}
\details{
 Compare methods facilitate quality assessment of synthetic data by comapring 
 them with the original observed data sets. The data themselves (for class 
 \code{synds}) or models fitted to them (for class \code{fit.synds}) are 
 compared.   
}
\value{
 The value returned by \code{compare} depends on the class of its argument. 
 See the documentation of the particular methods for details. 
}
\seealso{
  \code{\link{compare.synds}}, \code{\link{compare.fit.synds}}
}
