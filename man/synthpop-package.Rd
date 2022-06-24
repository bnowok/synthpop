\name{synthpop-package}
\alias{synthpop-package}
\alias{synthpop}
\docType{package}
\title{Generating synthetic versions of sensitive microdata for statistical disclosure control}
\description{Generate synthetic versions of a data set using parametric or CART methods.}
\details{
\tabular{ll}{
Package: \tab synthpop\cr
Type: \tab Package\cr
Version: \tab 1.7-1\cr
Date: \tab 2022-06-24\cr
License: \tab GPL-2 | GPL-3\cr
}
Synthetic data are generated from the original (observed) data by the function 
\code{\link{syn}}. The package includes also tools to compare synthetic data with the 
observed data (\code{\link{compare.synds}}) and to fit (generalized) linear model to 
synthetic data (\code{\link{lm.synds}}, \code{\link{glm.synds}}) and compare the estimates 
with those for the observed data (\code{\link{compare.fit.synds}}). More extensive 
documentation with illustrative examples is provided in the package vignette.}
\author{
Beata Nowok, Gillian M Raab, and Chris Dibben based on package \pkg{mice} 
(2.18) by Stef van Buuren and Karin Groothuis-Oudshoorn

Maintainer: Beata Nowok <beata.nowok@gmail.com>
}
\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.
}
\keyword{package}
