\name{codebook.syn}
\alias{codebook.syn}
\title{Makes a codebook from a data frame}
\description{
  Describes features of variables in a data frame relevant for synthesis.
}
\usage{
codebook.syn(data, maxlevs = 3)
}
\arguments{
  \item{data}{a data frame with a data set to be synthesised.}
  \item{maxlevs}{the number of factor levels above which separate tables with 
    all labels are returned as part of \code{labs} component.}
}
\value{
  A list with two components.
  
  \code{tab} - a data frame with the following information about each variable:
  \item{name}{variable name}
  \item{class}{class of variable}
  \item{nmiss}{number of missing values (\code{NA})}
  \item{perctmiss}{percentage of missing values}
  \item{ndistinct}{number of distinct values (excluding missing values)}
  \item{details}{range for numeric variables, maximum length for character 
    variables, labels for factors with <= maxlevs levels}

 \code{labs} - a list of extra tables with labels for each factor with number 
   of levels greater than \code{maxlevs}.
}
\examples{
codebook.syn(SD2011)
}


