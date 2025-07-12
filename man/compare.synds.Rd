\name{compare.synds}
\alias{compare.synds}
\alias{compare.data.frame}
\alias{compare.list}
\alias{print.compare.synds}
\title{
 Compare univariate distributions of synthesised and observed data
}
\description{
 Compare synthesised data set with the original (observed) data set
 using percent frequency tables and histograms. When more than one
 synthetic data set has been generated (\code{object$m > 1}), by 
 default pooled synthetic data are used for comparison.

 This function can be also used with synthetic data NOT created by
 \code{syn()}, but then an additional parameter \code{cont.na} might
 need to be provided.
}
\usage{
\method{compare}{synds}(object, data, vars = NULL,
        msel = NULL, stat = "percents", breaks = 20, ngroups =5,
        nrow = 2, ncol = 2, rel.size.x = 1,
        utility.stats = c("pMSE", "S_pMSE", "df"),
        utility.for.plot = "S_pMSE",
        cols = c("#1A3C5A","#4187BF"),
        plot = TRUE, table = FALSE, 
        print.flag = TRUE, \dots)

\method{compare}{data.frame}(object, data, vars = NULL, cont.na = NULL,
        msel = NULL, stat = "percents", breaks = 20,ngroups = 5,
        nrow = 2, ncol = 2, rel.size.x = 1,
        utility.stats = c("pMSE", "S_pMSE", "df"),
        utility.for.plot = "S_pMSE",
        cols = c("#1A3C5A","#4187BF"),
        plot = TRUE, table = FALSE, 
        print.flag = TRUE, compare.synorig = TRUE, \dots)

\method{compare}{list}(object, data, vars = NULL, cont.na = NULL,
        msel = NULL, stat = "percents", breaks = 20,ngroups = 5,
        nrow = 2, ncol = 2, rel.size.x = 1,
        utility.stats = c("pMSE", "S_pMSE", "df"),
        utility.for.plot = "S_pMSE",
        cols = c("#1A3C5A","#4187BF"),
        plot = TRUE, table = FALSE, 
        print.flag = TRUE, compare.synorig = TRUE, \dots)

\method{print}{compare.synds}(x, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised
    data set'. It is typically created by function \code{syn()} and it includes
    \code{object$m} synthesised data set(s) as \code{object$syn}. Alternatively,
    when data are synthesised not using \code{syn()}, it can be a data frame
    with a synthetic data set or a list of data frames with synthetic data sets,
    all created from the same original data with the same variables and the same
    method.}
  \item{data}{an original (observed) data set.}
  \item{vars}{variables to be compared. If \code{vars} is \code{NULL}
    (the default) all synthesised variables are compared.}
  \item{cont.na}{a named list of codes for missing values for continuous
    variables if different from the \code{R} missing data code \code{NA}.
    The names of the list elements must correspond to the variables names for
    which the missing data codes need to be specified.}
  \item{msel}{index or indices of synthetic data copies for which a comparison
    is to be made. If \code{NULL} pooled synthetic data copies are compared
    with the original data.}
  \item{stat}{determines whether tables and plots present percentages
    \code{stat = "percents"}, the default, or counts \code{stat = "counts"}.
    If \code{m > 1} and \code{msel = NULL} average counts for synthetic data
    are derived.}
  \item{breaks}{the number of cells for the histogram.}
  \item{ngroups}{the number of groups used to categorise numeric 
  variables when calculating the one-way utility measures.}
  \item{nrow}{the number of rows for the plotting area.}
  \item{ncol}{the number of columns for the plotting area.}
  \item{rel.size.x}{a number representing the relative size of x-axis labels.}
  \item{utility.stats}{a single string or a vector of strings that determines
    which utility measures to print. Must be a selection from:  
    \code{"VW"}, \code{"FT"},\code{"JSD"}, \code{"SPECKS"}, \code{"WMabsDD"}, 
    \code{"U"}, \code{"G"}, \code{"pMSE"}, \code{"PO50"}, \code{"MabsDD"}, 
    \code{"dBhatt"},  \code{"S_VW"}, \code{"S_FT"}, \code{"S_JSD"}, 
    \code{"S_WMabsDD"}, \code{"S_G"}, \code{"S_pMSE"}, \code{"df"}. 
    If \code{utility.stats = "all"}, all of these will be printed. For more 
    information see the details section for \code{\link{utility.tab}}.}
  \item{utility.for.plot}{a single string that determines which utility 
    measure to print in facet labels of the plot. Set to \code{NULL} to
    print variable names only.}
  \item{cols}{bar colors.}
  \item{plot}{a logical value with default set to \code{TRUE} indicating
    whether plots should be produced.}
  \item{table}{a logical value with default set to \code{FALSE} indicating
    whether tables should be printed.}
\item{print.flag}{a logical value with default set to \code{TRUE} indicating
    whether a message should be printed as the statistic for each 
    variable are calculated.}
  \item{compare.synorig}{a logical value  to determine if the functions
  \code{synorig.compare()} should be used to check that data sets can be 
   compared. Used when the synthetic data are supplied as a data.frame or 
   a list when  default set to TRUE.} 
  \item{\dots}{additional parameters.}
  \item{x}{an object of class \code{compare.synds}.}
}
\details{
Missing data categories for numeric variables are plotted on the same plot
as non-missing values. They are indicated by \code{miss.} suffix.

Numeric variables with fewer than 6 distinct values are changed to factors
in order to make plots more readable.
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
\seealso{
  \code{\link{multi.compare}}
}
\examples{
ods <- SD2011[ , c("sex", "age", "edu", "marital", "ls", "income")]
s1  <- syn(ods, cont.na = list(income = -8))

### synthetic data provided as a 'synds' object
compare(s1, ods, vars = "ls")
compare(s1, ods, vars = "income", stat = "counts",
        table = TRUE, breaks = 10)

### synthetic data provided as 'data.frame'
compare(s1$syn, ods, vars = "ls")
compare(s1$syn, ods, vars = "income", cont.na = list(income = -8),
        stat = "counts", table = TRUE, breaks = 10)
}
