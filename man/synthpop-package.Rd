\name{synthpop-package}
\alias{synthpop-package}
\alias{synthpop}
\docType{package}
\title{Generating synthetic versions of sensitive microdata for statistical disclosure control}
\description{Generate high- or low-fidelity synthetic versions of a data set using conditional models
   based on parametric or CART methods, or methods that fit the joint distributions of categorical 
   variables (catall and ipf).}
\details{
\tabular{ll}{
Package: \tab synthpop\cr
Type: \tab Package\cr
Version: \tab 1.9-3\cr
Date: \tab 2025-11-05\cr
License: \tab GPL-2 | GPL-3\cr
}
Synthetic data are generated from the original (observed) data by the function 
\code{\link{syn}}. The package includes also tools to compare synthetic data with the 
observed data (\code{\link{compare.synds}}) and to fit (generalized) linear model to 
synthetic data (\code{\link{lm.synds}}, \code{\link{glm.synds}}) and compare the estimates 
with those for the observed data (\code{\link{compare.fit.synds}}). More extensive 
documentation on how to create synthetic data, with illustrative examples, is provided
in the package vignette \bold{synthpop}. Since that vignette was written more
methods have been added to synthpop, including mthods for categorical variables
based on log-linear models that can be made differentially private.

Now the package also includes functions to eavaluate the utility and
disclosure risk of synthetic data. For details see the vignettes
\bold{utility} and \bold{disclosure}. You can access all the vignettes via the index link at the bottom of 
this help page (\code{\link{synthpop-package}})
}

\author{
Beata Nowok, Gillian M Raab, and Chris Dibben 
}

\references{
  Elliot, M. (2014) Final report on the disclosure risk associated with the synthetic data produced by the SYLLS team. Report 2015-2, Cathie Marsh Centre for Census and Survey Research (CCSR).
  
  Nowok, B. Utility of synthetic microdata generated using tree-based methods (2015) Paper presented at the Privacy in Statistical Databases Conference 2016; Dubrovnik, Croatia, 14-16 September 2016 .
    
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke   creation of synthetic data in R. \bold{Journal of Statistical Software}, \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.
  
  Raab, G.M., Nowok, B., and Dibben, C. (2016) Practical data synthesis for large samples \bold{Journal of Privacy and Confidentiality}, \bold{7}(3):67-97. \doi{10.29012/jpc.v7i3.407}. 

  Raab, G.M., Nowok, B., and Dibben, C. (2016) Guidelines for producing useful synthetic data  \doi{10.48550/arXiv.1712.04078} An earlier version was presented at the Privacy in Statistical Databases Conference 2016; Dubrovnik, Croatia, 14-16 September 2016 

  Nowok, B., Raab, G.M. and Dibben, C. (2017) Providing bespoke synthetic data for the UK Longitudinal Studies and other sensitive data with the synthpop package for R \bold{Statistical Journal of the IAOS, 33(3):785-796}. \doi{10.3233/SJI-150153}.

  Raab, G.M., Nowok, B., and Dibben, C. (2021) Assessing, visualizing and improving the utility of synthetic data. Available at\doi{10.48550/arXiv.2109.12717}. An earlier version was presented at the Joint UNECE/Eurostat expert meeting on statistical data confidentiality; Poznan, Poland, 1-3 December 2021.
  
  Raab, G.M. (2022) Utility and Disclosure Risk for Differentially Private Synthetic Categorical Data, Chapter in Privacy in Statistical Databases 2022. Published in Springer Series Lecture notes in Computer Science.  Also available at \doi{10.48550/arXiv.2206.01362}.  
  
  Raab, G.M., Nowok, B., and Dibben, C. (2024) Practical privacy metrics for synthetic data, Vignette in synthpop package. Also available at \doi{10.48550/arXiv.2406.16826}. 
   
  Raab, G.M. (2024) Privacy risk from synthetic data: practical proposals. Chapter in Privacy in Statistical databases 2024. published in Springer Series Lecture notes in Computer Science.  Also available at  \doi{10.48550/arXiv.2409.04257}. 

}
\keyword{package}
