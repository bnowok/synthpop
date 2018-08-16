\name{compare.synds}
\alias{compare.synds}
\alias{print.compare.synds}
\title{
 Compare univariate distributions of synthesised and observed data
}
\description{
 Compare synthesised data set with the original (observed) data set 
 using percent frequency tables and histograms. When more than one 
 synthetic data set has been generated (\code{object$m>1}), by default 
 pooled synthetic data are used for comparison.
}
\usage{
\method{compare}{synds}(object, data, vars = NULL, msel = NULL, 
  breaks = 20, nrow = 2, ncol = 2, rel.size.x = 1, 
  cols = c("#1A3C5A","#4187BF"), stat = "percents", ...)

\method{print}{compare.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{syn()} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{data}{an original (observed) data set.}
  \item{vars}{variables to be compared. If \code{vars} is \code{NULL} 
    (the default) all synthesised variables are compared.}
  \item{msel}{index or indices of synthetic data copies for which a comparison 
    is to be made. If \code{NULL} pooled synthetic data copies are compared 
    with the original data.}    
  \item{breaks}{the number of cells for the histogram.}  
  \item{nrow}{the number of rows for the plotting area.}
  \item{ncol}{the number of columns for the plotting area.}
  \item{rel.size.x}{a number representing the relative size of x-axis labels.}
  \item{cols}{bar colors.}
  \item{stat}{determines whether tables and plots present percentages 
   \code{stat = "percents"}, the default, or counts \code{stat = "counts"}. 
   If \code{m > 1} and \code{msel = NULL} average counts for synthetic data 
   are presented.}
  \item{\dots}{additional parameters.}

  \item{x}{an object of class \code{compare.synds}.}
}
\details{
Missing data categories for numeric variables are plotted on the same plot 
as non-missing values. They are indicated by \code{miss.} suffix.
}
\value{
An object of class \code{compare.synds} which is a list including a list 
of comparative frequency tables (\code{tables}) and a ggplot object
(\code{plots}) with bar charts/histograms. If multiple plots are produced
they and their corresponding frequency tables are stored as a list. 
}
\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.
}
\examples{
ods <- SD2011[ , c("sex","age","edu","marital","ls","income")]
s1  <- syn(ods)
compare(s1, ods, vars = "ls")
compare(s1, ods, vars = "income", stat = "counts", breaks = 10)
}
