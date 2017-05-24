\name{utility.tab}
\alias{utility.tab}
\alias{print.utility.tab}
\title{Tabular utility}
\description{
  Produce tables from observed and synthesized data and calculates 
  utility measures to compare them with their expectation if the synthesising 
  model is correct.}
\usage{
utility.tab(object, data, vars = NULL, ngroups = 5, with.missings = TRUE,
            digits = 2, print.tables = FALSE, print.zdiff = FALSE, \dots) 
       
\method{print}{utility.tab}(x, digits = x$digits, print.tables = x$print.tables, print.zdiff = x$print.zdiff, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} or 
    \code{syn.strata()} and it includes \code{object$m} number of synthesised 
    data set(s), as well as \code{object$syn} the synthesised data set, 
    if \code{m = 1}, or a list of \code{m} such data sets.}
  \item{data}{the original (observed) data set.}
  \item{vars}{a single string or a vector of strings with the names of 
    variables to be used to form the table.}
  \item{ngroups}{if numerical (non-factor) variables are included they will be 
    classified into this number of groups to form tables. For data values that 
    are all unique the groups are formed from the percentiles of the 
    distribution. But when there are repeated values identical observations are 
    always in the same group. This may produce a smaller number of groups than 
    specified with \code{ngroups} in some cases.}
  \item{with.missings}{a logical value indicating whether the table is to 
    include rows and columns for missing data for categorical variables.}
  \item{\dots}{additional parameters.}   
  \item{x}{an object of class \code{utility.tab}.}
  \item{digits}{an integer indicating the number of decimal places 
    for printing statistics and \code{tab.Zscore}.}
  \item{print.tables}{a logical value that determines if tables of observed, 
    synthesised are to be printed.} 
  \item{print.zdiff}{a logical value that determines if tables of  Z scores for 
  differences between observed and expected are to be printed.} 
}
\details{Forms tables of observed and synthesised values for the variables 
  specified in \code{vars}. Two utility measures are calculated from the cells 
  of the tables, a measure of fit proposed by Voas and Williams 
  \code{sum((observed-synthesied)^2/[(observed + synthesised)/2)])} and one 
  proposed by Freeman and Tukey \code{4*sum((observed^(0.5)-synthesised^(0.5))^2))}. 
  In both cases those cells where observed and synthesised are both zero do not 
  contribute to the sum. If the synthesising model is correct both of these 
  measures should have chi-square distributions for large samples.
}
\value{An object of class \code{utility.tab} which is a list with the following 
  components: 
  \item{m}{Number of synthetic data sets in object, i.e. \code{object$m}.}
  \item{tab.obs}{a table from the observed data.} 
  \item{UtabFT}{a vector with \code{object$m} values for the Freeman Tukey 
    utility measure.} 
  \item{UtabVW}{a vector with \code{object$m} values for the Voas Williamson 
    utility measure.} 
  \item{df}{a vector of degrees of freedom for the chi-square tests which equal 
    to one minus the number of cells in the table with any observed or 
    synthesised counts.} 
  \item{ratioFT}{a vector with ratios of \code{UtabFT} to \code{df}.}
  \item{ratioVW}{a vector with ratios of \code{UtabVW} to \code{df}.} 
  \item{stdFT}{a vector with standardised \code{UtabFT} in standard deviation 
    units compared to its null distribution if the synthesising model is correct.}
  \item{stdVW}{a vector with standardised \code{UtabVW} in standard deviation 
    units compared to its null distribution if the synthesising model is correct.}
  \item{nempty}{a vector of length \code{object$m} with number of cells 
    not contributing to the statistics.} 
  \item{with.missings}{whether missing values were included in the table.}
  \item{nsyn.missings}{number of missing observations that would contribute to 
    the tables from each synthetic data sets.}
  \item{nobs.missings}{number of missing observations that would contribute to 
    the tables from the observed data.}
  \item{tab.obs}{a table from the observed data.} 
  \item{tab.syn}{a table or a list of \code{m} tables from the synthetic data.}  
  \item{tab.zdiff}{a table or a list of \code{m} tables of Z statistics for 
    differences between observed and synthesised cells of the tables. Large 
    absolute values indicate a large contribution to lack-of-fit.}    
}

\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.

  Read, T.R.C. and Cressie, N.A.C. (1988) \emph{Goodness--of--Fit Statistics for 
  Discrete Multivariate Data}, Springer--Verlag, New York.  
  
  Voas, D. and Williamson, P. (2001) Evaluating goodness-of-fit measures for 
  synthetic microdata. \emph{Geographical and Environmental Modelling},
  \bold{5}(2), 177-200.
}

\seealso{
 \code{\link{utility.gen}}
}

\examples{
ods <- SD2011[1:1000, c("sex", "age", "edu", "marital")]

s1 <- syn(ods, m = 10)
utility.tab(s1, ods, vars = c("marital", "sex"))

s2 <- syn(ods, m = 1)
u2 <- utility.tab(s2, ods, vars = c("marital", "age"), groups = TRUE, 
                  ngroups = 3, with.missings = FALSE)
print(u2, print.zdiff = TRUE)
}
