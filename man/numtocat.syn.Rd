\name{numtocat.syn}
\alias{numtocat.syn}
\title{Group numeric variables before synthesis}
\description{
Selected numeric variables are grouped into factors with ranges selected from 
the data.
}
\usage{
numtocat.syn(data, numtocat = NULL, print.flag = TRUE, cont.na = NULL, 
             catgroups = 5, style.groups = "quantile") 
}
\arguments{
  \item{data}{a data frame.}
  \item{numtocat}{a vector of numbers or variable names of numeric variables 
    to be grouped into factors. If \code{NULL} all the numeric variables in data 
    will be grouped.}
  \item{print.flag}{if TRUE a list of grouped variables is printed.}
  \item{cont.na}{a named list that gives the values of the named variables to be 
    treated as separate categories, often missing values like \code{-8}. See the 
    corresponding parameter of \code{syn()}.}
  \item{catgroups}{a single integer or a vector of integers indicating the target 
    number of groups for the variables in numtocat in the same order as numtocat, 
    or as their relative postions in data. The achieved number of groups may be 
    different if, for example there are fewer than \code{ngroups} distinct values.}
  \item{style.groups}{parameter of the function \code{classInt()} that determines 
    how the breaks used to categorise each variable are chosen. See the help file 
    for \code{classInt()} for details. The default setting \code{"quantile"} makes 
    groups of approximately equal size. To divide into approximately equal ranges 
    we suggest using \code{"fisher"}.}
  }
\value{
A list with the following components:
  \item{data}{a data frame with the numeric variables replaced by factors grouped 
    into ranges.}
  \item{breaks}{a named list of the breaks used to divide each numeric variable 
    into categories.}
  \item{levels}{a named list of the levels for the categories of each numeric 
    variable.}
  \item{orig}{a data frame with the original numeric data.}
  \item{cont.na}{a named list of the levels for the categorical version of each 
    numeric variable.}
  \item{numtocat}{names of the variables changed to categories.}
  \item{ind}{positions in data of  the variables changed to categories.}
  }
\examples{
SD2011.cat <- numtocat.syn(SD2011, cont.na = list(income = -8 , unempdur = -8, 
nofriend = -8))
summary(SD2011.cat$data)
}  

  



