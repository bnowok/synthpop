\name{utility.tab}
\alias{utility.tab}
\alias{utility.tab.synds}
\alias{utility.tab.data.frame}
\alias{utility.tab.list}
\alias{print.utility.tab}
\title{Tabular utility}
\description{
  Produces tables from observed and synthesised data and calculates 
  utility measures to compare them with their expectation if the 
  synthesising model is correct.
  
  It can be also used with synthetic data NOT created by \code{syn()}, 
  but then an additional parameter\code{cont.na} might need to be provided.
}
\usage{
\method{utility.tab}{synds}(object, data, vars = NULL, ngroups = 5, 
                            useNA = TRUE, print.tables = length(vars) < 4, 
                            print.stats = "VW", print.zdiff = FALSE, 
                            digits = 3, k.syn = FALSE, \dots) 

\method{utility.tab}{data.frame}(object, data, vars = NULL, cont.na = NULL, 
                                 ngroups = 5, useNA = TRUE, 
                                 print.tables = length(vars) < 4, 
                                 print.stats = "VW", print.zdiff = FALSE, 
                                 digits = 3, k.syn = FALSE, \dots) 
                            
\method{utility.tab}{list}(object, data, vars = NULL, cont.na = NULL, 
                           ngroups = 5, useNA = TRUE, 
                           print.tables = length(vars) < 4, 
                           print.stats = "VW", print.zdiff = FALSE, 
                           digits = 3, k.syn = FALSE, \dots) 

                            
\method{print}{utility.tab}(x, print.tables = NULL, 
                            print.zdiff = NULL, print.stats = NULL, 
                            digits = NULL, \dots)
}

\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} or 
    \code{syn.strata()} and it includes \code{object$m} number of synthesised 
    data set(s), as well as \code{object$syn} the synthesised data set, 
    if \code{m = 1}, or a list of \code{m} such data sets. Alternatively, 
    when data are synthesised not using \code{syn()}, it can be a data frame 
    with a synthetic data set or a list of data frames with synthetic data sets, 
    all created from the same original data with the same variables and the same 
    method.}
  \item{data}{the original (observed) data set.}
  \item{vars}{a single string or a vector of strings with the names of 
    variables to be used to form the table.}
  \item{cont.na}{a named list of codes for missing values for continuous
    variables if different from the \code{R} missing data code \code{NA}.
    The names of the list elements must correspond to the variables names for 
    which the missing data codes need to be specified.}    
  \item{ngroups}{if numerical (non-factor) variables are included they will be 
    classified into this number of groups to form tables. Classification is 
    performed using \code{classIntervals()} function for \code{n = ngroups}. 
    By default, \code{style = "quantile"} to get appropriate groups for skewed 
    data. Problems for variables with a small number of unique values are handled 
    by selecting only unique values of breaks. Arguments of \code{classIntervals()} 
    may be, however, specified in the call to \code{utility.tab()}.}
  \item{useNA}{determines if NA values are to be included in tables.}
  \item{print.tables}{a logical value that determines if tables of observed and 
    synthesised are to be printed.} 
  \item{print.stats}{a single string or a vector of strings that determines 
    which utility measures to print to compare the observed and synthetic tables: 
    "VW" for Voas Williams, "FT" for Freeman Tukey, "G" for a modified likelihood 
    ratio for tables, "JSD" for Jensen-Shannon distance, and "pct.correct" for
    percent correctly predicted. If \code{print.stats = "all"}, all of the measures
    mentioned above will be printed. For more information see the details section 
    below.}   
  \item{print.zdiff}{a logical value that determines if tables of  Z scores for 
   differences between observed and expected are to be printed.} 
  \item{digits}{an integer indicating the number of decimal places 
    for printing statistics, \code{tab.zdiff} and mean results for \code{m > 1}.
    Note that for p-value it is fixed at 4 and cannot be changed.}
  \item{k.syn}{a logical indicator as to whether the sample size itself has 
    been synthesised. The default value is \code{FALSE}, which will apply 
    to synthetic data created by synthpop.}    
  \item{\dots}{additional parameters; can be passed to classIntervals() function.}   
  \item{x}{an object of class \code{utility.tab}.}
}
\details{Forms tables of observed and synthesised values for the variables 
  specified in \code{vars}. Several utility measures are calculated from the cells 
  of the tables, as described below. The one we recommend is a measure of fit 
  proposed by Voas and Williams \code{sum((observed-synthesied)^2/[(observed + synthesised)/2)])}.
  Others are one proposed by Freeman and Tukey \code{4*sum((observed^(0.5)-synthesised^(0.5))^2))} 
  and a modified likelihood ratio for tables (\code{G}). This last statistic is not 
  recommended for the case where there will be any zero entries in either the real 
  or observed data. In all cases those cells where observed and synthesised are both 
  zero do not contribute to the sum. If the synthesising model is correct all of these 
  measures should have chi-square distributions for large samples.
  Two other measures calculated are the Jensen-Shannon distance (\code{JSD}) and 
  the percentage of the observations in the combined tables (\code{pct.correct}) 
  where the majority of observations in each grouping are in agreement with category 
  (real or synthetic) of the observation.
}
\value{An object of class \code{utility.tab} which is a list with the following 
  components: 
  \item{m}{number of synthetic data sets in object, i.e. \code{object$m}.}
  \item{UtabFT}{a vector with \code{object$m} values for the Freeman Tukey 
    utility measure.} 
  \item{UtabVW}{a vector with \code{object$m} values for the Voas Williamson 
    utility measure.} 
  \item{UtabG}{a vector with \code{object$m} values for the adjusted likelihood ratio 
    utility measure.} 
  \item{df}{a vector of degrees of freedom for the chi-square tests which equal 
    to the number of cells in the tables with any observed or 
    synthesised counts minus one when \code{k.syn == FALSE} or equal to the 
    the number of cells when \code{k.syn == TRUE}.} 
  \item{ratioFT}{a vector with ratios of \code{UtabFT} to \code{df}.}
  \item{pvalFT}{a vector with \code{object$m} p-values for the chi-square 
    tests for the Freeman Tukey utility measure.}
  \item{ratioVW}{a vector with ratios of \code{UtabVW} to \code{df}.} 
  \item{pvalVW}{a vector with \code{object$m} p-values for the chi-square 
    tests for the Voas Williamson utility measure.}
  \item{ratioG}{a vector with ratios of \code{UtabVW} to \code{df}.} 
  \item{pvalG}{a vector with \code{object$m} p-values for the chi-square 
    tests for the adjusted likelihood ratio utility measure.}
  \item{JSD}{a vector with \code{object$m} values for the Jensen-Shannon distance.} 
  \item{pct.correct}{a vector with \code{object$m} values for percent correctly 
    predicted.} 
  \item{nempty}{a vector of length \code{object$m} with number of cells 
    not contributing to the statistics.} 
  \item{tab.obs}{a table from the observed data.} 
  \item{tab.syn}{a table or a list of \code{m} tables from the synthetic data.}  
  \item{tab.zdiff}{a table or a list of \code{m} tables of Z statistics for 
    differences between observed and synthesised cells of the tables. Large 
    absolute values indicate a large contribution to lack-of-fit.}
  \item{digits}{an integer indicating the number of decimal places 
    for printing statistics, \code{tab.zdiff} and mean results for \code{m > 1}.}
  \item{print.tables}{a logical value that determines if tables of observed and 
    synthesised are to be printed.} 
  \item{print.stats}{a single string or a vector of strings with utility measures 
   to be printed out.}   
  \item{print.zdiff}{a logical value that determines if tables of  Z scores for 
   differences between observed and expected are to be printed.} 
  \item{n}{number of observation in the original dataset.}
  \item{k.syn}{a logical indicator as to whether the sample size itself has 
    been synthesised.}
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
ods <- SD2011[1:1000, c("sex", "age", "marital", "nofriend")]

s1 <- syn(ods, m = 10, cont.na = list(nofriend = -8))
utility.tab(s1, ods, vars = c("marital", "sex"), print.stats = "all")

s2 <- syn(ods, m = 1, cont.na = list(nofriend = -8))
u2 <- utility.tab(s2, ods, vars = c("marital", "age", "sex"), ngroups = 3)
print(u2, print.tables = TRUE, print.zdiff = TRUE)

### synthetic data provided as 'data.frame'
utility.tab(s2$syn, ods, vars = c("marital", "nofriend"), ngroups = 3, 
            print.tables = TRUE, cont.na = list(nofriend = -8), digits = 4)
}
