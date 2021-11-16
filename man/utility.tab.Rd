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
  but then an additional parameter \code{cont.na} might need to be provided.
}
\usage{
\method{utility.tab}{synds}(object, data, vars = NULL, ngroups = 5,
            useNA = TRUE, max.table = 1e6,
            print.tables = length(vars) < 4,
            print.stats = c("pMSE", "S_pMSE", "df"),
            print.zdiff = FALSE, print.flag = TRUE,
            digits = 4, k.syn = FALSE, \dots)

\method{utility.tab}{data.frame}(object, data, vars = NULL, cont.na = NULL,
            ngroups = 5, useNA = TRUE, max.table = 1e6,
            print.tables = length(vars) < 4,
            print.stats = c("pMSE", "S_pMSE", "df"),
            print.zdiff = FALSE, print.flag = TRUE,
            digits = 4, k.syn = FALSE, \dots)

\method{utility.tab}{list}(object, data, vars = NULL, cont.na = NULL,
            ngroups = 5, useNA = TRUE, max.table = 1e6,
            print.tables = length(vars) < 4,
            print.stats = c("pMSE", "S_pMSE", "df"),
            print.zdiff = FALSE, print.flag = TRUE,
            digits = 4, k.syn = FALSE, \dots)


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
  \item{max.table}{a maximum table size. You could try increasing the default 
    value, but memory problems are likely.}
  \item{ngroups}{if numerical (non-factor) variables are included they will be
    classified into this number of groups to form tables. Classification is
    performed using \code{classIntervals()} function for \code{n = ngroups}.
    By default, \code{style = "quantile"} to get appropriate groups for skewed
    data. Problems for variables with a small number of unique values are handled
    by selecting only unique values of breaks. Arguments of \code{classIntervals()}
    may be, however, specified in the call to \code{utility.tab()}.}
  \item{useNA}{determines if NA values are to be included in tables.}
  \item{print.tables}{a logical value that determines if tables of observed and
    synthesised data are to be printed. By default tables are printed if they 
    have up to three dimensions.}
  \item{print.stats}{a single string or a vector of strings that determines
    which utility measures to print. Must be a selection from:  
     \code{"VW"}, \code{"FT"},\code{"JSD"}, \code{"SPECKS"}, \code{"WMabsDD"}, 
     \code{"U"}, \code{"G"}, \code{"pMSE"}, \code{"PO50"}, \code{"MabsDD"}, 
     \code{"dBhatt"},  \code{"S_VW"}, \code{"S_FT"}, \code{"S_JSD"}, 
     \code{"S_WMabsDD"}, \code{"S_G"}, \code{"S_pMSE"}, \code{"df"}, \code{dfG}. 
    If \code{print.stats = "all"}, all of these will be printed. For more 
    information see the details section below.}
  \item{print.zdiff}{a logical value that determines if tables of Z scores for
   differences between observed and expected are to be printed.}
  \item{print.flag}{a logical value that determines if messages are to be printed 
    during computation.}
  \item{digits}{an integer indicating the number of decimal places for printing
    statistics, \code{tab.zdiff} and mean results for \code{m > 1}.}
  \item{k.syn}{a logical indicator as to whether the sample size itself has
    been synthesised. The default value is \code{FALSE}, which will apply
    to synthetic data created by synthpop.}
  \item{\dots}{additional parameters; can be passed to classIntervals() function.}
  \item{x}{an object of class \code{utility.tab}.}
}
\details{Forms tables of observed and synthesised values for the variables
specified in \code{vars}. Several utility measures are calculated from the cells
of the tables, as described below. Details of all of these measures can be found 
in Raab et al. (2021). If the synthesising model is correct the measures
\code{VW}, \code{FT}, \code{G} and \code{JSD} should have chi-square distributions 
with \code{df} degrees of freedom for large samples. Standardised versions of each 
measure are available (e.g. \code{S_VW} for \code{VW}, where \code{S_VW = VW/df}) 
that will have an expected value of \code{1} if the synthesising model is correct. 
Four other measures are calculated by considering the table as a prediction model. 
The propensity score mean-squared error \code{pMSE}, and from a comparison of 
propensity scores for the synthetic and original data the Kolmogorov-Smirnov 
statistic \code{SPECKS} and the Wilcoxon rank-sum statistic \code{U} and also 
the percentage of the observations correctly predicted in the combined tables over 
50\%(\code{PO50}) where the majority of observations in each grouping are in 
agreement with category (real or synthetic) of the observation. The first of these 
\code{pMSE} is identical except for a constant to \code{VW}. No expected values are 
computed for the last three of these measures, but they can be obtained by replication 
from \code{utility.gen()}. 
Three further measures are calulated from the tables. The mean absolute difference 
in distributions: firstly \code{MabsDD}, the avarage absolute difference in the 
proportions of original and synthetic data from all the cells in the table. 
Secondly a weighted version of this measure \code{WMabsDD} where the weights are 
proportional to the inverse of the variance of the absolute differences so that 
this measure can be standardised by its expected value, \code{df}. Finally the 
Bhattacharyya distances \code{BhattD} derived from the overlap of the histograms 
of the original and synthetic data sets. 
}
\value{An object of class \code{utility.tab} which is a list with the following
  components:
  \item{m}{number of synthetic data sets in object, i.e. \code{object$m}.}
  \item{VW}{a vector with \code{object$m} values for the Voas Williamson
    utility measure.; linearly related to \code{pMSE}.}
  \item{FT}{a vector with \code{object$m} values for the Freeman-Tukey
    utility measure.}
  \item{JSD}{a vector with \code{object$m} values for the Jensen-Shannaon 
    divergence for comparing the tables.}
  \item{SPECKS}{a vector with \code{object$m} values for the Kolmogorov-Smirnov 
    statistic for comparing the propensity scores for the original and synthetic 
    data.}
  \item{WMabsDD}{a vector with \code{object$m} values of the weighted mean 
    absolute difference in distributions for original and synthetic data.}
  \item{U}{a vector with \code{object$m} values of the Wilcoxon statistic 
    comparing the propensity scores for the original and synthetic data.}
  \item{G}{a vector with \code{object$m} values for the adjusted likelihood 
    ratio utility measure.}
  \item{pMSE}{a vector with \code{object$m} values of the propensity score 
    mean-squared error; linearly related to \code{VW}.} 
  \item{PO50}{a vector with \code{object$m} values of the percentage over 
    50\% of observations correctly predicted from the propensity scores 
    linearly related to \code{SPECKS} and \code{MabsDD}.} 
  \item{MabsDD}{a vector with \code{object$m} values of the mean absolute 
    difference in distributions for original and synthetic data linearly 
    related to \code{SPECKS} and \code{PO50}.}
  \item{dBhatt}{a vector with \code{object$m} values of the Bhattacharyya 
    distances between the synthetic and original data, linearly related to 
    the square root of \code{FT}.}
  \item{S_VW}{\code{VW/df}.}
  \item{S_FT}{\code{FT/df}.}
  \item{S_JSD}{ \code{JSD}/df.}
  \item{S_WMabsDD}{WMabsDD/df.}
  \item{S_G}{\code{G/df}.}
  \item{S_pMSE}{standardised measure from \code{pMSE}, identical to \code{S_VW}.}
  \item{df}{a vector of degrees of freedom for the chi-square tests which equal
    to the number of cells in the tables with any observed or
    synthesised counts minus one when \code{k.syn == FALSE} or equal to the
    the number of cells when \code{k.syn == TRUE}.}
  \item{dfG}{degrees of freedom used in standardising \code{G}.}  
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

  Raab, G.M., Nowok, B. and Dibben, C. (2021). Assessing, visualizing and 
  improving the utility of synthetic data. Available from 
  \url{https://arxiv.org/abs/2109.12717}.

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
