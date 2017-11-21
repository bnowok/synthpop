\name{syn.satcat}
\alias{syn.satcat}
\title{Synthesis from a saturated model based on all combinations of the 
  predictor variables.}
\description{Synthesises one variable (\code{y}) from all possible 
  combinations of its precitors (\code{x}). A bootstrap sample is created 
  from the original values of \code{y} within each unique combinations of 
  of \code{xp} (the syntheisied values of the grouping variable).}
\usage{
syn.satcat(y, x, xp, proper = FALSE, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n} for the satcat variable.}
  \item{x}{a matrix (\code{n} x \code{p}) with the original predictor variables
    for \code{y}.}
  \item{xp}{a matrix (\code{k} x \code{p}) with synthetic values of \code{x}.}
  \item{proper}{if \code{proper = TRUE} values are a bootstrapped 
    sample from the synthesised data.} 
  \item{\dots}{additional parameters.}
}
\details{It is intended that the variables in \code{x} are categorical (factor) 
variables. If  \code{y} is also a categorical variable \code{syn.satcat} will 
give the same results as fitting a saturated polychotomous regression model but 
will usually be much faster. \code{syn.satcat} will fail with an error message 
if previous syntheses have generated a combination of variables in \code{xp} 
that was not present in \code{x}.}
\value{A vector of length \code{k} with synthetic values of \code{y}.}
\examples{
ods <- SD2011[, c("region", "sex", "agegr", "placesize")]

s1 <- syn(ods, method = c("sample", "cart", "satcat", "cart"))

### mostly fails because too many small categories
s2 <- syn(ods, method = c("sample", "cart", "cart", "satcat"))  
}
\keyword{datagen}
