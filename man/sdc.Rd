\name{sdc}
\alias{sdc}
\title{Tools for statistical disclosure control (sdc)}
\description{Labeling, top and bottom coding, smoothing numeric data, and
removing different types of unique records  defined by keys from synthetic data.
The function calls \code{\link{replicated.uniques}} to identify the rows
to be excluded from the synthetic data set(s)}
\usage{
sdc(object, data,keys = NULL, prefix = NULL, suffix = NULL, label = NULL, 
rm.uniques.in.orig = FALSE, rm.replicated.uniques = FALSE, 
recode.vars = NULL, bottom.top.coding = NULL, 
 recode.exclude = NULL, smooth.vars = NULL)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{the original (observed) data set.}
  \item{keys}{Variables to be used as quasi-identifiers to check for unique
  combinations. Passed to  \code{\link{replicated.uniques}} to exclude rows
  in the synthetic data.}
  \item{prefix}{A character string to be added as a prefix to all variable 
    names in the synthetic data set(s)}
  \item{suffix}{A character string to be added as a suffix to all variable 
    names in the synthetic data set(s)}
  \item{label}{a single string with a label to be added to the synthetic data 
    sets as a new variable to make it clear that the data are synthetic/fake.}
  \item{rm.uniques.in.orig}{a logical value indicating whether unique 
    replicates of key variables that are present in the orginal data set 
    should  be removed from synthetic data set(s).}
  \item{rm.replicated.uniques}{a logical value indicating whether unique 
    replicates of key variables that are also unique in the orginal data set 
    should  be removed.}
  \item{recode.vars}{a single string or a vector of strings with name(s) 
    of variable(s) to be bottom- or/and top-coded.}
  \item{bottom.top.coding}{a list of two-element vectors specifing 
    bottom and top codes for each variable in \code{recode.vars}. 
    If there is no need for one of bottom or top coding \code{NA} should be used.
    If only one variable is to be recoded, codes can be given as 
    a two-element vector.}
  \item{recode.exclude}{a list specifying for each variable in 
    \code{recode.vars} values to be excluded from recoding, e.g. missing data 
    codes. If all non-missing values should be considered for recoding \code{NA} should be 
    used if missing values are present. If only one variable is to be recoded, code(s) can be given as 
    a single number or a vector.}
  \item{smooth.vars}{a single string or a vector of strings with name(s) 
    of numeric variable(s) to be smoothed (\code{\link{smooth.spline}} function is used).}
}
\value{
  An \code{object} provided as an argument adjusted in accordance with the 
  other parameters' values.   
}
\seealso{
  \code{\link{replicated.uniques}}
}
\examples{
ods <- SD2011[1:1000,c("sex","age","region","edu","marital","income")]
s1 <- syn(ods, m = 2)
s1.sdc <- sdc(s1, ods, keys = c("sex","age","region"),suffix = "_synthetic",
label="false_data", rm.uniques.in.orig = TRUE,
recode.vars = c("age","income"),
bottom.top.coding = list(c(20,80),c(NA,2000)),
recode.exclude = list(NA,c(NA,-8)))
head(s1.sdc$syn[[2]])
}
\keyword{manip}