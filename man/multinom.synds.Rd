\name{multinom.synds}
\alias{multinom.synds}
\title{Fitting multinomial models to synthetic data}
\description{
  Fits multinomial models to the synthesised data set(s) 
  using the \code{\link[nnet]{multinom}} function.
}
\usage{
multinom.synds(formula, data,  ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be estimated. 
    A typical model has the form \code{response ~ predictors}. See
    the documentation of \code{\link[nnet]{multinom}} and
    \code{\link[stats]{formula}} for details.}
  \item{data}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{\link{syn}} and it includes \code{data$m} 
    synthesised data set(s).}
  \item{\dots}{additional parameters passed to \code{\link[nnet]{multinom}}.}
}
\value{
  To print the results the print function (\code{\link{print.fit.synds}}) can
  be used. The \code{\link{summary}} function (\code{\link{summary.fit.synds}}) 
  can be used to obtain the combined results of models fitted to each of the 
  \code{m} synthetic data sets.
  
  An object of class \code{fit.synds}. It is a list with the following 
  components:
  \item{call}{the original call to \code{multinom.synds}.}
  \item{mcoefavg}{combined (average) coefficient estimates.}
  \item{mvaravg}{combined (average) variance estimates of \code{mcoef}.}
  \item{analyses}{an object summarising the fit to each synthetic data set
    or a list of \code{m} such objects. Note that this is different from 
    the object created by \code{summary.multinom} to make it compatible 
    with other fitting methods. In particular the coefficients are vectors, 
    not matrices.}
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
  \code{\link[nnet]{multinom}}, \code{\link[synthpop]{glm.synds}}, 
  \code{\link[synthpop]{polr.synds}}, \code{\link[synthpop]{print.fit.synds}},
  \code{\link[synthpop]{summary.fit.synds}}, \code{\link[synthpop]{compare.fit.synds}}
}
\examples{
ods <- SD2011[1:1000, c("sex", "age", "edu", "marital", "ls", "smoke")]
s1 <- syn(ods, m = 3)
f1 <- multinom.synds(edu  ~ sex + age, data = s1)
summary(f1)
print(f1, msel = 1:2)
compare(f1, SD2011)
}
\keyword{multivariate}

