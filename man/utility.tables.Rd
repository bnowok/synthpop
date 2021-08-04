\name{utility.tables}
\alias{utility.tables}
\alias{utility.tables.synds}
\alias{utility.tables.data.frame}
\alias{utility.tables.list}
\alias{print.utility.tables}
\title{Tables and plots of utility measures}
\description{
  Calculates and plots tables of utility measures. The calculations of 
  utility measures are done by the functions \code{\link{utility.tab}}
  or \code{\link{utility.gen}}. Options are all one-way tables, 
  all two-way tables or three-way tables for a specified third variable
  along with pairs of all other variables.

  This function can be also used with synthetic data NOT created by 
  \code{syn()}, but then an additional parameters \code{not.synthesised} 
  and \code{cont.na} might need to be provided.  
}

\usage{
\method{utility.tables}{synds}(object, data,
                           tables = "twoway", method = "tab", 
                           vars = NULL, third.var = NULL,
                           useNA = TRUE, ngroups = 6,
                           utility.stats = NULL, plot.stat = NULL,
                           convergence.action = "setNA", plot = TRUE,
                           print.tab.res = FALSE, digits.tab.res = 2,
                           maxit = 100, max.scale = NULL, plottitle = NULL,
                           nworst = 5, ntabstoprint  = 0, k.syn = FALSE, \dots)
               
\method{utility.tables}{data.frame}(object, data, 
                           cont.na = NULL, not.synthesised = NULL, 
                           tables = "twoway", method = "tab", 
                           vars = NULL, third.var = NULL, 
                           useNA = TRUE, ngroups = 6, 
                           utility.stats = NULL, plot.stat = NULL,
                           convergence.action = "setNA", plot = TRUE,   
                           print.tab.res = FALSE, digits.tab.res = 2,
                           maxit = 100, max.scale = NULL, plottitle = NULL,  
                           nworst = 5, ntabstoprint  = 0, k.syn = FALSE, \dots)          

\method{utility.tables}{list}(object, data, 
                           cont.na = NULL, not.synthesised = NULL, 
                           tables = "twoway", method = "tab", 
                           vars = NULL, third.var = NULL, 
                           useNA = TRUE, ngroups = 6, 
                           utility.stats = NULL, plot.stat = NULL,
                           convergence.action = "setNA", plot = TRUE,   
                           print.tab.res = FALSE, digits.tab.res = 2,
                           maxit = 100, max.scale = NULL, plottitle = NULL,  
                           nworst = 5, ntabstoprint  = 0, k.syn = FALSE, \dots)

\method{print}{utility.tables}(x, print.tab.res = NULL, digits.tab.res = NULL, 
                               plot = NULL, plottitle = NULL, max.scale = NULL, 
                               nworst = NULL, ntabstoprint = NULL, \dots)}

\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s) as \code{object$syn}. This a single 
    data set when \code{object$m = 1} or a list of length \code{object$m} when 
    \code{object$m > 1}. Alternatively, when data are synthesised not using 
    \code{syn()}, it can be a data frame with a synthetic data set or a list 
    of data frames with synthetic data sets, all created from the same original 
    data with the same variables and the same method.}  
  \item{data}{the original (observed) data set.}
  \item{cont.na}{a named list of codes for missing values for continuous
    variables if different from the \code{R} missing data code \code{NA}.
    The names of the list elements must correspond to the variables names for 
    which the missing data codes need to be specified.} 
  \item{not.synthesised}{a vector of variable names for any variables that has 
    been left unchanged in the synthetic data.}
  \item{tables}{defines the type of tables to produce. For method \code{"gen"} 
    it must be "oneway" or "twoway". For method \code{"tab"} options are 
    \code{"oneway"}, \code{"twoway"} (default) or \code{"threeway"}. 
    If set to \code{"oneway"} or \code{"twoway"} all possible tables from 
    \code{vars} are produced. For \code{"threeway"}, \code{third.var} must be 
    specified and all three-way tables between this variable and other pairs of 
    variables are produced.}
  \item{method}{this determines which method is used to compute utility metrics. 
    It can be \code{"tab"} for tabular utility or \code{"gen"} for general 
    utility.}
  \item{vars}{a vector of strings with the names of variables to be used to form 
    the table, or a vector of variable numbers in the original data. Defaults to 
    all variables in both original and synthetic data.}    
  \item{third.var}{required when \code{tables} is \code{"threeway"} a variable 
    to make the third variable with all other pairs}    
  \item{useNA}{determines if \code{NA} values are to be included in tables. Only 
    applies for method \code{"tab"}.}
  \item{ngroups}{if numerical (non-factor) variables included with 
    \code{method = "tab"} will be classified into this number of groups to form
    tables. Classification is performed using \code{classIntervals()} function 
    for \code{n = ngroups}. By default, \code{style = "quantile"}, to get 
    appropriate groups for skewed data. Problems for variables with a small 
    number of unique  values are handled by selecting only unique values of 
    breaks. Arguments of \code{classIntervals()} may be, however, specified 
    in the call to \code{utility.tables()}.}
  \item{utility.stats}{statistics to include in the table of results. Choices for 
    \code{method = "tab"} are any combination of \code{"VW"}, \code{"FT"}, 
    \code{"G"}, \code{"JSD"}, \code{"pct.correct"}, or \code{"df"}. 
    For \code{method = "gen"} they are \code{"pMSE"}, \code{"pMSEExp"}, 
    \code{"utilR"}, \code{"pval"}, \code{"df"}, or \code{"pct.correct"}. 
    See \code{\link{utility.tab}} and \code{\link{utility.gen}} for 
    explanations of measures.}
  \item{plot.stat}{statistics to plot. Choice for \code{method = "tab"} is 
    \code{"VW"}, \code{"FT"}, \code{"G"}, \code{"JSD"} or \code{"pct.correct"}.  
    For \code{method = "gen"} choice is \code{"utilR"} or \code{"pct.correct"}.
    If \code{"pcorrect"} is chosen the statistic plotted will be the \% correctly
    predicted in excess of 50\%, since 50\% represents no ability to predict 
    when an observation is synthetic. See \code{\link{utility.tab}} and 
    \code{\link{utility.gen}} for explanations of measures.}
  \item{convergence.action}{for \code{method = "gen"} determines what happens 
    when a model fails to converge in the number of iterations specified. 
    Possible values are \code{"none"} no action or \code{"setNA"} when the 
    values in the outputs for these tables and plots are set to \code{NA}.}
  \item{plot}{determines if plot will be produced when the result is printed.} 
  \item{print.tab.res}{logical value that determines if table of results is 
    to be printed.}  
  \item{digits.tab.res}{number of digits to print for table, except for 
    p-values that are always printed to 4 places.}  
  \item{maxit}{the maximum number of iterations allowed before a model for 
    \code{method} is considered not to have converged.}
  \item{max.scale}{a numeric value for the maximum value used in calculating 
    the shading of the plots. If it is \code{NULL} then the maximum value  
    will be replaced by the maximum value in the data.} 
  \item{plottitle}{title for the plot.} 
  \item{nworst}{...}
  \item{ntabstoprint}{...}
  \item{k.syn}{a logical indicator as to whether the sample size itself has 
    been synthesised.}
  \item{\dots}{additional parameters}   
  \item{x}{an object of class \code{utility.tables}.}
}

\details{Calculates tables of observed and synthesised values for the variables 
  specified in \code{vars} with the function \code{\link{utility.tab}} or
  \code{\link{utility.gen}} and produces tables and plots of one-way, two-way or 
  three-way utility tables formed from \code{vars}. Several options for utility 
  measures can be selected for printing or plotting. Details are in help files 
  for \code{\link{utility.tab}} and \code{\link{utility.gen}}. 
  
  The tables and variables with the worst utility scores are identified. 
  Visualisations of the matrices of utility scores are plotted. For threeway
  tables a third variable can be defined to select all tables involving that 
  variable for plotting. If it is not specified the variable with tables 
  giving the worst utility is selected as the third variable.
}

\value{An object of class \code{utility.tab} which is a list with the following 
  components: 
  \item{tab.res}{a table with all the selected measures for all combinations of 
    variables defined by \code{tables}, \code{third.var}, and \code{vars}.}
  \item{method}{see above.}  
  \item{plot.stat}{measure used in \code{mat} and \code{toplot}.} 
  \item{tables}{see above.} 
  \item{third.var}{see above.} 
  \item{toplot}{table used to produce the plot when printing.}
  \item{any.not.converged}{TRUE/FALSE to indicate if any general utility models 
    have not converged}
  \item{convergence.action}{see above.}
  \item{var.scores}{...} 
  \item{print.tab.res}{see above.}  
  \item{digits.tab.res}{see above.} 
  \item{plottitle}{see above.} 
  \item{max.scale}{see above.}
  \item{ntabstoprint}{...}  
  \item{nworst}{...}  
  \item{worstn}{...}   
  \item{worsttabs}{...}
}

\references{
  Read, T.R.C. and Cressie, N.A.C. (1988) \emph{Goodness--of--Fit Statistics 
  for Discrete Multivariate Data}, Springer--Verlag, New York.  
  
  Voas, D. and Williamson, P. (2001) Evaluating goodness-of-fit measures for 
  synthetic microdata. \emph{Geographical and Environmental Modelling},
  \bold{5}(2), 177-200.
}

\seealso{
 \code{\link{utility.tab}}, \code{\link{utility.gen}}
}

\examples{
ods <- SD2011[1:1000, c("sex", "age", "edu", "marital", "region", "income")]
s1 <- syn(ods)

### synthetic data provided as a 'synds' object  
(t1 <- utility.tables(s1, ods, utility.stats = "all", print.tab.res = TRUE))
### synthetic data provided as a 'data.frame' object
(t1 <- utility.tables(s1$syn, ods, utility.stats = "all", print.tab.res = TRUE))
t1$toplot

t2 <- utility.tables(s1, ods, tables = "oneway")
print(t2, max.scale = 3)

(t3 <- utility.tables(s1, ods, utility.stats = "all", tables = "threeway", 
                      third.var = "sex", print.tab.res = TRUE))

(t4 <- utility.tables(s1, ods, utility.stats = "all", tables = "threeway", 
                      third.var = "sex", useNA = FALSE, print.tab.res = TRUE))

(t5 <- utility.tables(s1, ods, method = "gen", utility.stat = "all", 
                      print.tab.res = TRUE))
t5$toplot

t6 <- utility.tables(s1, ods, method = "gen", tables = "oneway")
print(t6, max.scale = 3)
}
