\name{sdc}
\alias{sdc}
\title{Tools for statistical disclosure control (sdc)}
\description{Labeling and removing unique replicates of unique actual (observed) individuals.}
\usage{
sdc(object, data, label = NULL, rm.replicated.uniques = FALSE, 
 uniques.exclude = NULL, recode.vars = NULL, bottom.top.coding = NULL, 
 recode.exclude = NULL, smooth.vars = NULL)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{the original (observed) data set.}
  \item{label}{a single string with a label to be added to the synthetic data 
    sets as a new variable to make it clear that the data are synthetic/fake.}
  \item{rm.replicated.uniques}{a logical value indicating whether unique 
    replicates of units that are unique also in the orginal data set should 
    be removed.}
  \item{uniques.exclude}{a single string or a vector of strings with name(s) 
    of variable(s) to be excluded from the identification of uniques.}
  \item{recode.vars}{a single string or a vector of strings with name(s) 
    of variable(s) to be bottom- or/and top-coded.}
  \item{bottom.top.coding}{a list of two-element vectors specifing 
    bottom and top codes for each variable in \code{recode.vars}. 
    If there is no need for bottom or top coding \code{NA} should be used.
    If only one variable is to be recoded, codes can be given as 
    a two-element vector.}
  \item{recode.exclude}{a list specifying for each variable in 
    \code{recode.vars} values to be excluded from recoding, e.g. missing data 
    codes. If all values should be considered for recoding \code{NA} should be 
    used. If only one variable is to be recoded, code(s) can be given as 
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
ods <- SD2011[1:1000,c("sex","age","edu","marital","income")]
s1 <- syn(ods, m = 2)
s1.sdc <- sdc(s1, ods, label="false_data", rm.replicated.uniques = TRUE,
recode.vars = c("age","income"),
bottom.top.coding = list(c(20,80),c(NA,2000)),
recode.exclude = list(NA,c(NA,-8)))
}
\keyword{manip}