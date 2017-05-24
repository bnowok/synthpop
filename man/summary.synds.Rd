\name{summary.synds}
\alias{summary.synds}
\alias{print.summary.synds}
\title{Synthetic data object summaries}
\description{
Produces summaries of the synthesised variables. When more than one 
synthetic data set has been generated (object$m>1), by default summaries 
are calculated by averaging summary values for all synthetic data copies
(see \code{msel} argument).  
}
\usage{
\method{summary}{synds}(object, msel = NULL, maxsum = 7, 
digits = max(3, getOption("digits")-3), \dots)

\method{print}{summary.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}; a result of a call to 
    \code{\link{syn}}.}
  \item{msel}{index or indices of synthetic data copies for which a summary 
    is desired. If \code{NULL} (default) summaries are calculated by averaging 
    summary values for all synthetic data copies.}
  \item{maxsum}{integer, indicating how many levels should be shown for factors.}
  \item{digits}{integer, used for number formatting with \code{\link{format}}.}    
  \item{\dots}{additional arguments passed to \code{\link{summary}}.}
  \item{x}{an object of class \code{summary.synds}.}
}
\details{See \code{\link[base]{summary}} for more details.}
\value{
  An object of class \code{summary.synds}, which is a list with the following 
  components:
  \item{m}{the number of synthetic versions of the original (observed) data.}
  \item{msel}{index or indices of synthetic data copies for which a summary 
    is produced. If \code{NULL} summaries are calculated by averaging 
    summary values for all synthetic data copies.}
  \item{method}{a vector of synthesising methods applied to each variable
    in the saved synthesised data.}
  \item{result}{a table or a list of tabels (if more than one synthetic data 
    set is selected) with summaries of synthesised variables.}
}
\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.
}
\seealso{
\code{\link{summary}},\code{\link{print}}
}
\examples{
s1 <- syn(SD2011[,c("sex","age","edu","marital")], m = 3)
summary(s1)
summary(s1, msel = c(1,3))
}
