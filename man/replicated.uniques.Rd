\name{replicated.uniques}
\alias{replicated.uniques}
\title{Replications in synthetic data}
\description{
  Determines which unique units in the synthesised data set(s)
  replicates unique units in the original observed data set.
}
\usage{
replicated.uniques(object, data, exclude = NULL)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{the original observed data set.}
  \item{exclude}{a single string or a vector of strings with name(s) 
    of variable(s) to be excluded from the identification of uniques.}
}
\value{
  A list with the following components:
  \item{replications}{a vector (for \code{object$m = 1}) or a data frame
   with \code{object$m} columns (for \code{object$m > 1}) with logical 
   values indicateing duplicates in \code{m}th synthetic data set.}
  \item{no.replications}{a single number or a vector of \code{object$m} integers
   indicating the number of duplicates in the synthetic data set(s).}  
  \item{no.uniques}{a number of unique individuals in the original data set.}  
  \item{per.replications}{a single number or a vector of \code{object$m} numeric 
   values indicating the percentage of duplicates in the synthetic data set(s).} 
}
\seealso{
  \code{\link{sdc}}
}
\examples{
ods <- SD2011[1:1000,c("sex","age","edu","marital","smoke")]
s1 <- syn(ods, m = 2)
replicated.uniques(s1,ods)
}


