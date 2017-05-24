\name{read.obs}
\alias{read.obs}
\title{Importing original data sets form external files}
\description{Imports data data sets form external files into a data frame. 
Currently supported files include: sav (SPSS), dta (Stata), xpt (SAS),
csv (comma-separated file), tab (tab-delimited file) and 
txt (delimited text files). For SPSS, Stata and SAS it uses functions from 
the \code{foreign} package with some adjustments where necessary.}
\usage{
read.obs(file, convert.factors = TRUE, lab.factors = FALSE, 
export.lab = FALSE, ...)
}
\arguments{
  \item{file}{the name of the file (including extension) which the data are 
   to be read from.}
  \item{convert.factors}{a logical value indicating whether variables with 
  value labels in Stata and SPSS should be converted into R factors with 
  those levels.}
  \item{lab.factors}{a logical value indicating whether variables with 
  complete value labels but imported using their numeric codes 
  (\code{convert.factors = FALSE}) should be converted 
  from numeric to factor variables.}
  \item{export.lab}{a logical variable indicating whether labels from SPSS or 
  Stata should be exported to an external file.}
  \item{\dots}{additional parameters passed to read functions.}
}
\value{A data frame with an imported data set. For SPSS, Stata and SAS it has 
attributes with labels.}
\seealso{
 \code{\link{write.syn}}
}
