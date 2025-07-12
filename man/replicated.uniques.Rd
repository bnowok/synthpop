\name{replicated.uniques}
\alias{replicated.uniques}
\alias{print.repuniq.synds}
\title{Replications in synthetic data}
\description{
  Determines which unique units in the synthesised data set(s)
  have  combinations of variables in the keys as follows:
  
  1) unique in original data 
  
  2) unique in the synthetic data set(s)
  
  3) unique in synthetic data and present,but not necessarily unique in original
  
  4) unique in synthetic and unique in original.
  
  For each of 3) and 4) results are returned that identify the rows in the
  synthetic data with each type of unique.
  This function is called by \code{\link{sdc}} where there are options to
  include each type of unique.
}
\usage{
replicated.uniques(object, data, keys = names(data))

\method{print}{repuniq.synds}(x, \dots)
}


\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{the original observed data set.}
  \item{keys}{Variables to be used as quasi-identifiers to check for unique
  combinations.} 
  \item{\dots}{additional parameters} 
  \item{x}{an object of class \code{repuniq.synds}; a result of a call to \code{replicated.uniques()}.}
}
\value{
  A list of class "repuniq.synds" with the following components:
  \item{m}{number of synthetic data sets in object \code{object$m} }
  \item{n}{number of rows in data \code{object$n}}
  \item{k}{number of rows in of synthetic data set(s) in object \code{object$k}}
  \item{res_tab}{Table or list of tables with numbers and percentages of uniques}
  \item{synU.rm}{A vector of length \code{object$k} TRUE/FALSE values wherea
  TRUE value identifies a unique in synthetic and prtesent in the original} 
  \item{repU.rm}{A vector of length \code{object$k} TRUE/FALSE values where a
  TRUE value identifies a replicated unique} 
}
\seealso{
  \code{\link{sdc}}
}
\examples{
ods <- SD2011[1:1000,c("sex","age","region","edu","marital","smoke")]
s1 <- syn(ods, m = 2)
replicated.uniques(s1,ods, keys = c("sex","age","region"))
}


