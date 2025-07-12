\name{synorig.compare}
\alias{synorig.compare}
\alias{synorig.compare.list}
\title{
check synthetic and original if not produced by synthpop.
}
\description{
 Check, and attempt to adjust, synthetic datasets NOT created by  \code{syn()} 
 if not compatible with the original.
 The output is a list with 4 components, the first two giving adjusted
 versions of the input synthetic data and original and the third \code{needsfix}
 an indicator of whether the functions listed below are likely to
 run correctly on the adjusted data, the fourth \code{unchanged} indicates
 whether any changes have been made to the original or the synthetic data.
 
 Variables that are in both the synthetic and the original data are 
 checked to 1)  convert any character variables to R factors 2)  check
 that data types match 3)  check differences in whether variables
 have missing values 4)  check if the levels of factors agree and
 if not use the combination of both sets.
 
 \code{needsfix} becomes TRUE if 1) some variables are in the synthetic but
 not in the original 2) variables have different classes (after characters
 converted to factor) 3) there are missing values in the synthetic data
 but not in the original 4) some levels of a factor are in the synthetic data
 but not in the original.
Some warning messages are printed if level differences are not just due to
missing values.

The function can be run to compare original and synthetic or it can be called
by setting compare.synorig = TRUE when data are not supplied as an object
 of class synds created by synthpop in the following functions :
\code{utility.tab() utility.tables() compare() disclosure() disclosure.summary()}.
In this case the function attempts to correct differences and continue or, if this
is impossible, will prompt the user as to which variables need changing.
}
\usage{
synorig.compare(syn,orig, print.flag = TRUE)
}

\arguments{
  \item{syn}{A data set containing the synthesised data, or a list of such data
  sets. When \code{syn} is a list only the first member of the list is used and
  the syn component of the result is the adjusted version of this first member.}
    \item{orig}{The original data set. }
  \item{print.flag}{If TRUE prints non-essential summary messages.}
 }
\details{
Error messages explain briefly what adjustments have been made to the data
sets, what could not be fixed and what might need to be checked.
Both \code{orig} and \code{syn} are made
    into simple data frames for comparison (e.g. if tibbles or matrices)
}
\value{
A list with 3 components 
\item{syn}{adjusted version of \code{syn} or of the first element of
the list \code{syn} if it is a list.}
\item{orig}{adjusted version of \code{orig}.}
\item{needsfix}{TRUE/FALSE as to whether the outputs need to be fixed 
before utility and disclosure functions could be used on them. }
\item{unchanged}{TRUE/FALSE indicating if the outputs of the function 
are unchanged from the outputs. }
}
\references{to add}
\seealso{
  \code{\link{utility.gen} \link{utility.tab} \link{utility.tables}
  \link{compare.synds} \link{disclosure.synds}}
}
\examples{
library(synthpop)
orig <- SD2011[1:2000,]
pretendsyn <- SD2011[2001:5000, 1:5]
orig[,1] <- as.character(orig[,1])
codebook.syn(orig[,1:5])
newdata <- synorig.compare(pretendsyn, orig)
codebook.syn(newdata$orig[,1:5])
}
