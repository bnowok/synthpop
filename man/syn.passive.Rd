\name{syn.passive}
\alias{syn.passive}
\title{Passive synthesis}
\description{Derives a new variable according to a specified 
  function of synthesised data.}
\usage{syn.passive(data, func)}
\arguments{
  \item{data}{a data frame with synthesised data.} 
  \item{func}{a \code{formula} specifying transformations on data. 
    It is specified as a string starting with \code{~}.}
}
\details{
  Any function of the synthesised data can be specified. Note that several operators such as 
  \code{+}, \code{-}, \code{*} and \code{^} have different meanings in \code{formula} syntax. 
  Use the identity function \code{I()} if they should be interpreted as arithmetic operators, 
  e.g. \code{"~I(age^2)"}.
  Function \code{syn()} checks whether the passive assignment is correct in the original data 
  and fails with a warning if this is not true. The variables synthesised passively can be
  used to predict later variables in the synthesis except when they are numeric variables 
  with missing data. A warning is produced in this last case.
  }
\value{
A list with two components:
  \item{res}{a vector of length \code{k} including the result of applying 
    the \code{formula}.}
  \item{fit}{a name of the method used for synthesis (\code{"passive"}).}
}
\references{
  Van Buuren, S. and Groothuis-Oudshoorn, K. (2011).
  \code{mice}: Multivariate Imputation by Chained Equations
  in \code{R}. \emph{Journal of Statistical Software},
  \bold{45}(3), 1-67. \doi{10.18637/jss.v045.i03}
}
\author{Gillian Raab, 2021 based on Stef van Buuren, Karin Groothuis-Oudshoorn, 2000}

\seealso{
  \code{\link{syn}}
} 
\examples{
### the examples shows how inconsistencies in the SD2011 data are picked up 
### by syn.passive()
ods <- SD2011[, c("height", "weight", "bmi", "age", "agegr")]
ods$hsq <- ods$height^2
ods$sex <- SD2011$sex
meth <- c("cart", "cart", "~I(weight / height^2 * 10000)",  
          "cart", "~I(cut(age, c(15, 24, 34, 44, 59, 64, 120)))", 
          "~I(height^2)", "logreg")

\dontrun{
### fails for bmi 
s1 <- syn(ods, method = meth, seed = 6756, models = TRUE)

### fails for agegr
ods$bmi <- ods$weight / ods$height^2 * 10000  
s2 <- syn(ods, method = meth, seed = 6756, models = TRUE) 

### fails because of wrong order 
ods$agegr <- cut(ods$age, c(15, 24, 34, 44, 59, 64, 120))
s3 <- syn(ods, method = meth, visit.sequence = 7:1, 
          seed = 6756, models = TRUE)  
}

### runs without errors
ods$bmi   <- ods$weight / ods$height^2 * 10000  
ods$agegr <- cut(ods$age, c(15, 24, 34, 44, 59, 64, 120))
s4 <- syn(ods, method = meth, seed = 6756, models = TRUE) 
### bmi and hsq do not predict sex because of missing values
s4$models$sex 

### hsq with no missing values used to predict sex  
ods2 <- ods[!is.na(ods$height),]
s5 <- syn(ods2, method = meth, seed = 6756, models = TRUE) 
s5$models$sex

### agegr with missing values used to predict sex because not numeric
ods3 <- ods
ods3$age[1:4] <- NA
ods3$agegr <- cut(ods3$age, c(15, 24, 34, 44, 59, 64, 120))
s6 <- syn(ods3, method = meth, seed = 6756, models = TRUE) 
s6$models$sex  
}

\keyword{datagen}
