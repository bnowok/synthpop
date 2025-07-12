\name{syn.satcat}
\alias{syn.satcat}
\title{Synthesis from a saturated model based on all combinations of the 
  predictor variables.}
\description{Synthesises one variable (\code{y}) from all possible 
  combinations of its predictors (\code{x}). A bootstrap sample is created 
  from the original values of \code{y} within each unique combinations of 
  of \code{xp} (the synthesisied values of the grouping variable).
  Note that only combinations of predictor variable levels that appear in the 
  original data can be in the synthetic data. 
  The related method 
  (\code{syn.catall}) overcomes this by adding a small prior probability 
  to all zero cells in the cross tabulation from the original data that are not 
  structural zeros. But it has the limitation of requiring a complete 
  cross tabulation of all the variables.}
\usage{
syn.satcat(y, x, xp, proper = FALSE, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n} for the satcat variable.}
  \item{x}{a matrix (\code{n} x \code{p}) with the original predictor variables
    for \code{y}.}
  \item{xp}{a matrix (\code{k} x \code{p}) with synthetic values of \code{x}.}
  \item{proper}{if \code{proper = TRUE} \code{x} and \code{y} are replaced 
    with a bootstrap sample before synthesis, thus effectively sampling from 
    the posterior distribution of the model, given the data.} 
  \item{\dots}{additional parameters.}
}
\details{It is intended that the variables in \code{x} are categorical (factor) 
variables. If \code{y} is also a categorical variable \code{syn.satcat} will 
give the same results as fitting a saturated polychotomous regression model but 
will usually be much faster. \code{syn.satcat} will fail with an error message 
if previous syntheses have generated a combination of variables in \code{xp} 
that was not present in \code{x}. Use of the \code{syn.catall} method for 
grouped variables can overcome this.}
\value{
A list with two components:
  \item{res}{a data frame of dimension \code{k x p} containing the synthesised data.}
  \item{fit}{the cross-tabulation of the original predictor variables.}
}
\examples{
ods <- SD2011[, c("region", "sex", "agegr", "placesize")]

s1 <- syn(ods, method =  "satcat", seed = 7856)
s2 <- syn(ods, method = c("sample", "cart", "satcat", "cart"), seed = 7856)

\dontrun{
### mostly fails because previous synthesis has produced 
### combinations not found in the original data
s3 <- syn(ods, method = c("sample", "cart", "cart", "satcat"), seed = 7856)}
}
\keyword{datagen}
