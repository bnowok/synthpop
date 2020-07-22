\name{glm.synds, lm.synds}
\alias{glm.synds}
\alias{lm.synds}
\alias{print.fit.synds}
\title{Fitting (generalized) linear models to synthetic data}
\description{
  Fits generalized linear models or simple linear models to the synthesised 
  data set(s) using \code{\link[stats]{glm}} and \code{\link[stats]{lm}} 
  function respectively.
}
\usage{
glm.synds(formula, family = "binomial", data,  ...)
lm.synds(formula, data, ...)

\method{print}{fit.synds}(x, msel = NULL, \dots) 
}
\arguments{
  \item{formula}{a symbolic description of the model to be estimated. 
    A typical model has the form \code{response ~ predictors}. See
    the documentation of \code{\link[stats]{glm}} and
    \code{\link[stats]{formula}} for details.}
  \item{family}{a description of the error distribution 
    and link function to be used in the model. See the documentation of 
    \code{\link[stats]{glm}} and \code{\link[stats]{family}} for details.}
  \item{data}{an object of class \code{synds}, which stands
    for 'synthesised data set'. It is typically created by
    function \code{\link{syn}} and it includes \code{data$m} 
    synthesised data set(s).}
  \item{\dots}{additional parameters passed to \code{\link[stats]{glm}} 
    or \code{\link[stats]{lm}}.}
  \item{x}{an object of class \code{fit.synds}.}    
  \item{msel}{index or indices of synthetic data copies for which coefficient
    estimates are to be displayed. If \code{NULL} (default) the combined 
    (average) coefficient estimates are printed.}  
}
\value{
  The \code{\link{summary}} function (\code{\link{summary.fit.synds}}) can be used 
  to obtain the combined results of models fitted to each of the \code{m} synthetic  
  data sets.

  An object of class \code{fit.synds}. It is a list with the following 
  components:
  \item{call}{the original call to \code{glm.synds} or \code{lm.synds}.}
  \item{mcoefavg}{combined (average) coefficient estimates.}
  \item{mvaravg}{combined (average) variance estimates of \code{mcoef}.}
  \item{analyses}{\code{summary.glm} or \code{summary.lm} object respectively
    or a list of \code{m} such objects.}
  \item{fitting.function}{function used to fit the model.}    
  \item{n}{a number of cases in the original data.}  
  \item{k}{a number of cases in the synthesised data.}  
  \item{proper}{a logical value indicating whether synthetic data were 
    generated using proper synthesis.}
  \item{m}{the number of synthetic versions of the observed data.}
  \item{method}{a vector of synthesising methods applied to each variable
    in the saved synthesised data.}
  \item{incomplete}{a logical value indicating whether any of the variables in 
    the model were not synthesised.} 
  \item{mcoef}{a matrix of coefficients estimates from all \code{m} syntheses.}
  \item{mvar}{a matrix of variance estimates from all \code{m} syntheses.}
}
\seealso{
  \code{\link[stats]{glm}}, \code{\link[stats]{lm}},  
  \code{\link[synthpop]{multinom.synds}}, \code{\link[synthpop]{polr.synds}}, 
  \code{\link[synthpop]{compare.fit.synds}}, \code{\link[synthpop]{summary.fit.synds}}
}
\examples{
### Logit model 
ods <- SD2011[1:1000, c("sex", "age", "edu", "marital", "ls", "smoke")]
s1 <- syn(ods, m = 3)
f1 <- glm.synds(smoke ~ sex + age + edu + marital + ls, data = s1, family = "binomial")
f1
print(f1, msel = 1:2)
  
### Linear model
ods <- SD2011[1:1000,c("sex", "age", "income", "marital", "depress")]
ods$income[ods$income == -8] <- NA
s2 <- syn(ods, m = 3)
f2 <- lm.synds(depress ~ sex + age + log(income) + marital, data = s2)
f2
print(f2,1:3)
}
\keyword{multivariate}

