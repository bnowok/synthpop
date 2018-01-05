\name{utility.tab}
\alias{utility.tab}
\alias{print.utility.tab}
\title{Tabular utility}
\description{
  Produce tables from observed and synthesized data and calculates 
  utility measures to compare them with their expectation if the synthesising 
  model is correct.}
\usage{
utility.tab(object, data, vars = NULL, ngroups = 5, 
            print.tables = TRUE, print.zdiff = FALSE, 
            digits = 2, \dots) 
       
\method{print}{utility.tab}(x, print.tables = x$print.tables, 
  print.zdiff = x$print.zdiff, digits = x$digits, \dots)
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
    classified into this number of groups to form tables. Classification is 
    performed using \code{classIntervals()} function for \code{n = ngroups}. 
    By default, to avoid problems for variables with a small number of unique 
    values, \code{style = "fisher"}. Arguments of \code{classIntervals()} may be, 
    however, specified in the call to \code{utility.tab()}.}
  \item{digits}{an integer indicating the number of decimal places 
    for printing statistics, \code{tab.zdiff} and mean results for \code{m > 1}.}  
  \item{print.tables}{a logical value that determines if tables of observed and 
    synthesised are to be printed.} 
  \item{print.zdiff}{a logical value that determines if tables of  Z scores for 
   differences between observed and expected are to be printed.} 
  \item{\dots}{additional parameters; can be passed to classIntervals() function.}   
  \item{x}{an object of class \code{utility.tab}.}
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
  \item{m}{number of synthetic data sets in object, i.e. \code{object$m}.}
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
  \item{pvalFT}{a vector with \code{object$m} p-values for the chi-square 
    tests for the Freeman Tukey utility measure.}
  \item{pvalVW}{a vector with \code{object$m} p-values for the chi-square 
    tests for the Voas Williamson utility measure.}
  \item{nempty}{a vector of length \code{object$m} with number of cells 
    not contributing to the statistics.} 
  \item{tab.obs}{a table from the observed data.} 
  \item{tab.syn}{a table or a list of \code{m} tables from the synthetic data.}  
  \item{tab.zdiff}{a table or a list of \code{m} tables of Z statistics for 
    differences between observed and synthesised cells of the tables. Large 
    absolute values indicate a large contribution to lack-of-fit.}    
  \item{n}{number of observation in the original dataset.}
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
utility.tab(s2, ods, vars = c("marital", "age"), ngroups = 3, print.tables = TRUE)
u2 <- utility.tab(s2, ods, vars = c("marital", "age"), style = "pretty")
print(u2, print.tables = TRUE, print.zdiff = TRUE)
}
