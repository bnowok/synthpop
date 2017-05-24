\name{write.syn}
\alias{write.syn}
\title{Exporting synthetic data sets to external files}
\description{Exports synthetic data set(s) from synthesised data set 
(\code{synds}) object to external files of selected format. 
Currently supported file formats include: SPSS, Stata, SAS, csv, tab,
rda, RData and txt. For SPSS, Stata and SAS it uses functions from 
the \code{foreign} package with some adjustments where necessary.
Information about the synthesis is written into a separate text file.  

NOTE: Currently numeric codes and labels can be preserved correctly only 
for SPSS files imported into R using \code{\link{read.obs}} function.   
}
\usage{
write.syn(object, filename, 
filetype = c("SPSS", "Stata", "SAS", "csv", "tab", "rda", "RData", "txt"), 
convert.factors = "numeric", data.labels = NULL, save.complete = TRUE, 
extended.info = TRUE, ...)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{\link{syn}} and it includes \code{object$m} 
    synthesised data set(s).}
  \item{filename}{the name of the file (excluding extension) which the 
    synthetic data are to be written into. For multiple synthetic data sets
    it will be used as a prefix folowed respectively by \code{_1}, \code{_...},
     \code{_m}.}
  \item{filetype}{a desired format of the output files.}
  \item{convert.factors}{a single string indicating how to handle factors in 
    Stata output files. The default value is set to \code{"numeric"} in order 
    to preserve the numeric codes from the original data. 
    See \code{\link{write.dta}} for other possible values.}
  \item{data.labels}{a list with variable labels and value labels.}  
  \item{save.complete}{a logical value indicating whether a complete 
    'synthesised data set' (\code{synds}) object should be saved into a file
    (\code{synobject_filename.RData}).}
  \item{extended.info}{a logical value indicating whether extended information
    should be saved into an information file.}
  \item{\dots}{additional parameters passed to write functions.}  
}
\value{File(s) with synthesised data set(s) and a text file with information 
  about synthesis are produced. Optionally a complete synthesised data set
  object is saved into \code{synobject_filename.RData} file.}
\seealso{
  \code{\link{read.obs}}
}


