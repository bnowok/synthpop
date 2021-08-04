\name{multi.compare}
\alias{multi.compare}
\title{
  Multivariate comparison of synthesised and observed data
}
\description{
  Graphical comparisons of a variable (\code{var}) in the synthesised data set 
  with the original (observed) data set within subgroups defined by the 
  variables in a vector \code{by}. \code{var} can be a factor or a continuous 
  variable and the plots produced will depend on the class of \code{var}. 
  The variables in \code{by} will usually be factors or variables with only 
  a few values.  
}
\usage{
multi.compare(object, data, var = NULL, by = NULL, msel = NULL, 
  barplot.position = "fill", cont.type = "hist", y.hist = "count", 
  boxplot.point = TRUE, binwidth = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s).}
  \item{data}{an original (observed) data set.}
  \item{var}{variable to be compared between observed and synthetic data within 
    subgroups.}
  \item{by}{variables to be tabulated or cross-tabulated to form groups.}
  \item{barplot.position}{type of barplot. The default \code{"fill"} gives       
    a single bar with the proportions in each group while \code{"dodge"} gives 
    side-by-side bars with the numbers in each category.}
  \item{cont.type}{default \code{"hist"} gives histograms and \code{"boxplot"} 
    gives boxplots.}
  \item{y.hist}{defines y scale for histograms - \code{"count"} is default;
    \code{"density"} gives proportions.}
  \item{boxplot.point}{default (\code{TRUE}) adds individual points to boxplots.}
  \item{msel}{numbers of synthetic data sets to be used - must be numbers in 
    the range \code{1:object$m} - defaults to \code{1:object$m}}
  \item{binwidth}{sets width of a bin for histograms.}
  \item{\dots}{additional parameters that can be supplied to \code{\link{ggplot}}.}
}
\value{
Plots as specified above. A table of the numbers in the subgroups is printed 
to the R console.

Numeric variables with fewer than 6 distinct values are changed to factors 
in order to make plots more readable.
}
\seealso{
  \code{\link{compare.synds}}, \code{\link{compare.fit.synds}}
}
\examples{
### default synthesis of selected variables
vars <- c("sex", "age", "edu", "smoke")
ods  <- na.omit(SD2011[1:1000, vars])
s1 <- syn(ods)

### categorical var
multi.compare(s1, ods, var = "smoke", by = c("sex","edu"))

### numeric var
multi.compare(s1, ods, var = "age", by = c("sex"), y.hist = "density", binwidth = 5)
multi.compare(s1, ods, var = "age", by = c("sex", "edu"), cont.type = "boxplot")
}

