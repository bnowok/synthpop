\name{polr.synds}
\alias{polr.synds}
\title{Fitting ordered logistic models to synthetic data}
\description{
  Fits  ordered logistic  models to the synthesised data set(s) 
  using the \code{\link[MASS]{polr}} function.
}
\usage{
polr.synds(formula, data, ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be estimated. A typical 
    model has the form \code{response ~ predictors}. See the documentation of 
    \code{\link[MASS]{polr}} and \code{\link[stats]{formula}} for details.}
  \item{data}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{\link{syn}} and it 
    includes \code{data$m} synthesised data set(s).}
  \item{\dots}{additional parameters passed to \code{\link[MASS]{polr}}.}
}
\value{
  To print the results the print function (\code{\link{print.fit.synds}}) can
  be used. The \code{\link{summary}} function (\code{\link{summary.fit.synds}})
  can be used to obtain the combined results of models fitted to each of the 
  \code{m} synthetic data sets.
  
  An object of class \code{fit.synds}. It is a list with the following 
  components:
  \item{call}{the original call to \code{polr.synds}.}
  \item{mcoefavg}{combined (average) coefficient estimates.}
  \item{mvaravg}{combined (average) variance estimates of \code{mcoef}.}
  \item{analyses}{an object summarising the fit to each synthetic data set
    or a list of \code{m} such objects. Note that this is different from  
    the object created by \code{summary.polr} to make it compatible with other 
    fitting methods for the \code{compare.fit.synds} and \code{summary.fit.synds} 
    functions. In particular the coefficients combine the items coefficients 
    and zeta from \code{summary.polr} to make a combined vector of coefficients.}
  \item{fitting.function}{function used to fit the model.}    
  \item{n}{a number of cases in the original data.}  
  \item{k}{a number of cases in the synthesised data.}  
  \item{proper}{a logical value indicating whether synthetic data were 
    generated using proper synthesis.}
  \item{m}{the number of synthetic versions of the observed data.}
  \item{method}{a vector of synthesising methods applied to each variable
    in the saved synthesised data.}
  \item{incomplete}{a logical value indicating whether the dependent variable
    in the model was not synthesised.} 
  \item{mcoef}{a matrix of coefficients estimates from all \code{m} syntheses.}
  \item{mvar}{a matrix of variance estimates from all \code{m} syntheses.}
}
\seealso{
  \code{\link[MASS]{polr}}, \code{\link[synthpop]{glm.synds}},  
  \code{\link[synthpop]{multinom.synds}}, \code{\link[synthpop]{print.fit.synds}},
  \code{\link[synthpop]{summary.fit.synds}}, \code{\link[synthpop]{compare.fit.synds}}
}
\examples{
ods <- SD2011[1:1000, c("sex", "age", "edu", "marital", "ls", "smoke")]
s1 <- syn(ods, m = 3)
f1 <- polr.synds(edu  ~ sex + age, data = s1)
summary(f1)
print(f1, msel = 1:2)
compare(f1, SD2011)
}
\keyword{multivariate}

