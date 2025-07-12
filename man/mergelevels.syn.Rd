\name{mergelevels.syn}
\alias{mergelevels.syn}
\title{Merge levels of factors in a data frame}
\description{
 merges levels of selected variables in a data frame either according to minimum numbers in a
 category or  according to user-defined rules.
}
\usage{
mergelevels.syn(data, vars = NULL, newlabel = FALSE, addNA = FALSE,
                    print.flag = FALSE, minsize = 10, merge.byhand =
                    FALSE, merge.details = NULL)

}

\arguments{
  \item{data}{An observed data set before synthesis.}
  \item{vars}{a vector of names or numbers for the variables for which categories are to be merged.
    defaults to all factors in data when sert to NULL}
  \item{minsize}{The minimum size that of combined categories when \code{merge.byhand} is FALSE}
  \item{newlabel}{When merge.byhand is FALSE when \code{newlabel = TRUE} causes the new factor
   levels to be labelled with a concatenation  of the labels of all the categories that have
    contributed.  \code{newlabel = FALSE}  labels the new category "Other"}
  \item{addNA}{Causes the NA category to be included when determining which groups are below
    \code{minsize} and also to allow NA values, designated by "NA",  to be combined with others
    in merge.details}
  \item{print.flag}{prints tables of variables before and after recoding.}
  \item{merge.byhand}{Uses the information in \code{merge.byhand}  to create groupings, ignoring 
   the parameters \code{minsize} and \code{vars}}.
   \item{merge.details}{A named list of variable names with names giving the names of the 
     variables that will have levels merged. Each item is a vector where the first item is
     the name for the new combined level and the other entries are the levels to be merged.
     If it exists already the levels with small counts will be added to it, otehrwise 
     a new level will be formed.}
  
  
  
}
\value{
 A data frame of the same size and structure as \code{data} with levels of selected variables merged.}

\seealso{
  \code{\link{sdc}}
}
\examples{
test <- SD2011[1:20]
data.mlevs1 <- mergelevels.syn(test, vars = c(3,5,18:20),minsize = 20, addNA = TRUE, 
              print.flag = TRUE, newlabel = TRUE)

mlevs <- list(agegr = c("60+", "60-64"  , "65+"),socprof = c("NEW","UNEMPLOYED","FARMER"))
data.mlevs2 <-mergelevels.syn(test, merge.byhand = TRUE, merge.details = mlevs, 
           addNA=TRUE, print.flag = TRUE)  


}


