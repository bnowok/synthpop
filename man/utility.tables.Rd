\name{utility.tables}
\alias{utility.tables}
\alias{utility.tables.synds}
\alias{utility.tables.data.frame}
\alias{utility.tables.list}
\alias{print.utility.tables}
\title{Tables and plots of utility measures}
\description{
  Calculates and plots tables of utility measures. The calculations of 
  utility measures are done by the function \code{\link{utility.tab}}.
  Options are all one-way tables, all two-way tables or three-way tables 
  for a specified third variable along with pairs of all other variables.

  This function can be also used with synthetic data NOT created by 
  \code{syn()}, but then an additional parameters \code{not.synthesised} 
  and \code{cont.na} might need to be provided.  
}

\usage{
\method{utility.tables}{synds}(object, data,
               tables = "twoway", maxtables = 5e4,
               vars = NULL, third.var = NULL,
               useNA = TRUE, ngroups = 5,
               tab.stats = c("pMSE", "S_pMSE", "df"), 
               plot.stat = "S_pMSE", plot = TRUE, max.table = 1e07,
               print.tabs = FALSE, digits.tabs = 4,
               max.scale = NULL, min.scale = 0, plot.title = NULL,
               nworst = 5, ntabstoprint = 0, k.syn = FALSE, 
               low = "grey92", high = "#E41A1C",
               n.breaks = NULL, breaks = NULL, print.flag = TRUE, \dots)
               
\method{utility.tables}{data.frame}(object, data, 
               cont.na = NULL, not.synthesised = NULL, 
               tables = "twoway", maxtables = 5e4,
               vars = NULL, third.var = NULL, 
               useNA = TRUE, ngroups = 5, 
               tab.stats = c("pMSE", "S_pMSE", "df"), 
               plot.stat = "S_pMSE", plot = TRUE, max.table = 1e07, 
               print.tabs = FALSE, digits.tabs = 4,
               max.scale = NULL, min.scale = 0, plot.title = NULL,  
               nworst = 5, ntabstoprint = 0, k.syn = FALSE,
               low = "grey92", high = "#E41A1C",
               n.breaks = NULL, breaks = NULL, 
               compare.synorig = TRUE, print.flag = TRUE,\dots)

\method{utility.tables}{list}(object, data,
               cont.na = NULL, not.synthesised = NULL, 
               tables = "twoway", maxtables = 5e4, 
               vars = NULL, third.var = NULL, 
               useNA = TRUE, ngroups = 5, 
               tab.stats = c("pMSE", "S_pMSE", "df"), 
               plot.stat = "S_pMSE", plot = TRUE,  max.table = 1e07,
               print.tabs = FALSE, digits.tabs = 4,
               max.scale = NULL, min.scale = 0, plot.title = NULL,  
               nworst = 5, ntabstoprint = 0, k.syn = FALSE,
               low = "grey92", high = "#E41A1C",
               n.breaks = NULL, breaks = NULL, 
               compare.synorig = TRUE, print.flag = TRUE,\dots)

\method{print}{utility.tables}(x, print.tabs = NULL, digits.tabs = NULL, 
      plot = NULL, plot.title = NULL, max.scale = NULL, min.scale = NULL, 
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
  \item{tables}{defines the type of tables to produce. Options are 
    \code{"oneway"}, \code{"twoway"} (default) or \code{"threeway"}. 
    If set to \code{"oneway"} or \code{"twoway"} all possible tables from 
    \code{vars} are produced. For \code{"threeway"}, \code{third.var} may be 
    specified and all three-way tables between this variable and other pairs of 
    variables are produced. If a third variable is not specified the function
    chooses the variable with the largest median utility measure for all three-way
    tables it contributes to.}
  \item{maxtables}{maximum number of tables that will be produced. If number of 
    tables is larger, then utility is only measured for a sample of size 
    \code{maxtables}. You cannot produce plots of twoway or three way tables from 
    sampled tables}.
  \item{vars}{a vector of strings with the names of variables to be used to form 
    the table, or a vector of variable numbers in the original data. Defaults to 
    all variables in both original and synthetic data.}    
  \item{third.var}{when \code{tables} is \code{"threeway"} a variable 
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
  \item{tab.stats}{statistics to include in the table of results. Must be 
    a selection from: \code{"VW"}, \code{"FT"},\code{"JSD"}, \code{"SPECKS"}, 
    \code{"WMabsDD"}, \code{"U"}, \code{"G"}, \code{"pMSE"}, \code{"PO50"}, 
    \code{"MabsDD"}, \code{"dBhatt"},  \code{"S_VW"}, \code{"S_FT"}, 
    \code{"S_JSD"}, \code{"S_WMabsDD"}, \code{"S_G"}, \code{"S_pMSE"}, 
    \code{"df"}, \code{dfG}. If \code{tab.stats = "all"}, all of these will 
    be included. See \code{\link{utility.tab}} for explanations of measures.}
  \item{plot.stat}{statistics to plot. Choice  is \code{"VW"}, \code{"FT"},
    \code{"JSD"}, \code{"SPECKS"},  \code{"WMabsDD"}, \code{"U"}, \code{"G"}, 
    \code{"pMSE"}, \code{"PO50"}, \code{"MabsDD"}, \code{"dBhatt"}, 
    \code{"S_VW"}, \code{"S_FT"}, \code{"S_JSD"}, \code{"S_WMabsDD"}, 
    \code{"S_G"}, \code{"S_pMSE"}. See \code{\link{utility.tab}} for 
    explanations of measures.}
  \item{plot}{determines if plot will be produced when the result is printed.} 
  \item{max.table}{Value of maximum number of cells allowed in a table by the function  \code{\link{utility.tab}} } 
  \item{print.tabs}{logical value that determines if table of results is 
    to be printed.}  
  \item{digits.tabs}{number of digits to print for table, except for 
    p-values that are always printed to 4 places.}  
  \item{max.scale}{a numeric value for the maximum value used in calculating 
    the shading of the plots. If it is \code{NULL} then the maximum value  
    will be replaced by the maximum value in the data.}
  \item{min.scale}{a numeric value for the minimum value used in calculating 
    the shading of the plots. If it is \code{NULL} then the minimum value  
    will be replaced by zero.}
  \item{plot.title}{title for the plot.} 
  \item{nworst}{a number of variable combinations with worst utility scores 
    to be printed.}
  \item{ntabstoprint}{a number of tables to print for observed and synthetic 
    data with the worst utility.}
  \item{k.syn}{a logical indicator as to whether the sample size itself has 
    been synthesised.}
  \item{low}{colour for low end of the gradient.}
  \item{high}{colour for high end of the gradient.} 
  \item{n.breaks}{a number of break points to create if breaks are not given 
    directly.}
  \item{breaks}{breaks for a two colour binned gradient.}
  \item{compare.synorig}{a logical value  to determine if the functions
  \code{synorig.compare()} should be used to check that data sets can be 
   compared. Used when the synthetic data are supplied as a data.frame or 
   a list when  default set to TRUE.} 
  \item{print.flag}{Allows printing of message as metrics are calculated 
  for each element of the table. Default is TRUE.}
  \item{\dots}{additional parameters}   
  \item{x}{an object of class \code{utility.tables}.}
}

\details{Calculates tables of observed and synthesised values for the variables 
  specified in \code{vars} with the function \code{\link{utility.tab}} and produces  
  tables and plots of one-way, two-way or 
  three-way utility measures formed from \code{vars}. Several options for utility 
  measures can be selected for printing or plotting. Details are in help file 
  for \code{\link{utility.tab}}. 
  
  The tables and variables with the worst utility scores are identified. 
  Visualisations of the matrices of utility scores are plotted. For threeway
  tables a third variable can be defined to select all tables involving that 
  variable for plotting. If it is not specified the variable with tables 
  giving the worst utility is selected as the third variable.
}

\value{An object of class \code{utility.tab} which is a list with the following 
  components: 
  \item{tabs}{a table with all the selected measures for all combinations of 
    variables defined by \code{tables}, \code{third.var}, and \code{vars}.}
  \item{plot.stat}{measure used in \code{mat} and \code{toplot}.} 
  \item{tables}{see above.} 
  \item{third.var}{see above.} 
  \item{utility.plot}{plot of the selected utility measure.}
  \item{var.scores}{an average of utility scores for all combinations with other
    variables.} 
  \item{plot}{see above.}
  \item{print.tabs}{see above.}  
  \item{digits.tabs}{see above.} 
  \item{plot.title}{see above.} 
  \item{max.scale}{see above.}
  \item{min.scale}{see above.}  
  \item{ntabstoprint}{see above.}  
  \item{nworst}{see above.}  
  \item{worstn}{variable combinations with \code{nworst} worst utility scores.}   
  \item{worsttabs}{observed and synthetic cross-tabulations for \code{worstn}.}
}

\references{
  Read, T.R.C. and Cressie, N.A.C. (1988) \emph{Goodness--of--Fit Statistics 
  for Discrete Multivariate Data}, Springer--Verlag, New York.  
  
  Voas, D. and Williamson, P. (2001) Evaluating goodness-of-fit measures for 
  synthetic microdata. \emph{Geographical and Environmental Modelling},
  \bold{5}(2), 177-200.
}

\seealso{
 \code{\link{utility.tab}}
}

\examples{
ods <- SD2011[1:1000, c("sex", "age", "edu", "marital", "region", "income")]
s1 <- syn(ods)

### synthetic data provided as a 'synds' object  
(t1 <- utility.tables(s1, ods, tab.stats = "all", print.tabs = TRUE))
### synthetic data provided as a 'data.frame' object
(t1 <- utility.tables(s1$syn, ods, tab.stats = "all", print.tabs = TRUE))

t2 <- utility.tables(s1, ods, tables = "twoway")
print(t2, max.scale = 3)

(t3 <- utility.tables(s1, ods, tab.stats = "all", tables = "threeway", 
                      third.var = "sex", print.tabs = TRUE))

(t4 <- utility.tables(s1, ods, tab.stats = "all", tables = "threeway", 
                      third.var = "sex", useNA = FALSE, print.tabs = TRUE))

(t5 <- utility.tables(s1, ods,  tab.stats = "all", 
                      print.tabs = TRUE))
}
