\name{compare.fit.synds}
\alias{compare.fit.synds}
\alias{print.compare.fit.synds}
\title{
 Compare model estimates based on synthesised and observed data
}
\description{
 The same model that was used for the synthesised data set is fitted to the 
 observed data set. The coefficients with confidence intervals for the 
 observed data is plotted together with their estimates from synthetic data. 
 When more than one synthetic data set has been generated (\code{object$m>1}) 
 combining rules are applied. Analysis-specific utility measures are used to
 evaluate differences between synthetic and observed data.
}
\usage{
\method{compare}{fit.synds}(object, data, plot = "Z", 
  print.coef = FALSE, return.plot = TRUE, plot.intercept = FALSE, 
  lwd = 1, lty = 1, lcol = c("#1A3C5A","#4187BF"), 
  dodge.height = .5, point.size = 2.5, incomplete = FALSE, 
  population.inference = FALSE, ci.level = 0.95, ...)

\method{print}{compare.fit.synds}(x, print.coef = x$print.coef, \dots)
}
\arguments{
  \item{object}{an object of type \code{fit.synds} created by fitting a model 
   to synthesised data set using function \code{\link{glm.synds}} or
   \code{\link{lm.synds}}.}
  \item{data}{an original observed data set.}
  \item{plot}{values to be plotted: \code{"Z"} (Z scores) or \code{"coef"}
   (coefficients).}
  \item{print.coef}{a logical value determining whether tables of 
   estimates for the original and synthetic data should be printed.}
  \item{return.plot}{a logical value indicating whether a confidence interval
   plot should be returned.}   
  \item{plot.intercept}{a logical value indicating whether estimates for 
   intercept should be plotted.}  
  \item{lwd}{the line type.}  
  \item{lty}{the line width.}
  \item{lcol}{line colours.}
  \item{dodge.height}{size of vertical shifts for confidence intervals to 
   prevent overlaping.} 
  \item{point.size}{size of plotting symbols used to plot point estimates of 
   coefficients.}
  \item{incomplete}{a logical value indicating whether the method of Reiter 
   (2003) for what he terms partially synthetic data should be used for 
   inference. It requires multiple syntheses with \code{m} greater than the 
   number of coefficients, ideally at least 5 more.} 
  \item{population.inference}{a logical value indicating whether intervals for 
   inference to population quantities, as decribed by Karr et al. (2006), 
   should be calculated and plotted. This option suppresses the lack-of-fit test 
   and the standardised differences since these are based on differences 
   standardised by the original interval widths.} 
  \item{ci.level}{Confidence interval coverage as a proportion.} 
  \item{\dots}{additional parameters passed to \code{\link{ggplot}}.}
  
  \item{x}{an object of class \code{compare.fit.synds}.}
}
\details{This function can be used to evaluate whether the method used for 
synthesis is appropriate for the fitted model. If this is the case the 
estimates from the  synthetic dataof what would be expected from the original 
data \code{xpct(Beta)} \code{xpct(Z)} should not differ from the estimates from 
the observed data (\code{Beta} and \code{Z}) by more than would be expected from 
the standard errors (\code{se(Beta)} and \code{se(Z)}). For more details see the 
vignette on inference.
}

\value{
 An object of class \code{compare.fit.synds} which is a list with the 
 following components:
  \item{call}{the original call to fit the model to the synthesised data set.}
  \item{coef.obs}{a data frame including estimates based on the observed
   data: coefficients (\code{Beta}), their standard errors (\code{se(Beta)})
   and Z scores (\code{Z}).}
  \item{coef.syn}{a data frame including (combined) estimates based on 
   the synthesised data: point estimates of observed data coefficients 
   (\code{B.syn}), standard errors of those estimates (\code{se(B.syn)}), 
   estimates of the observed standard errors (\code{se(Beta).syn}), Z scores 
   estimates (\code{Z.syn}) and their standard errors (\code{se(Z.syn)}). 
   Note that \code{se(B.syn)} and \code{se(Z.syn)} give the standard errors 
   of the mean of the \code{m} syntheses and can be made very small by 
   increasing \code{m} (see the vignette on inference for more details).}
 \item{coef.diff}{a data frame containing standardized differences between 
   the coefficients estimated from the original data and those calculated 
   from the combined synthetic data. The difference is standardized by dividing 
   by the estimated standard error of the fit from the original. The corresponding 
   p-values are calculated from a standard Normal distribution and represent the 
   probability of achieving differences as large as those found if the model use 
   for synthesis is compatible with the model that generated the original data.}  
 \item{mean.abs.std.diff}{Mean absolute standardized difference (over all 
   coefficients).} 
 \item{ci.overlap}{a data frame containing the percentage of overlap between 
   the estimated synthetic confidence intervals and the original sample 
   confidence intervals for each parameter. When \code{population.inference = TRUE} 
   overlaps are calculated as suggested by Karr et al. (2006). Otherwise 
   a simpler overlap measure with intervals of equal length is calculated.}
 \item{mean.ci.overlap}{Mean confidence interval overlap (over all 
   coefficients).}    
 \item{lack.of.fit}{lack-of-fit measure from all \code{m} synthetic data 
   sets combined, calculated as follows, when \code{incomplete.method = FALSE}. 
   The vector of mean differences (\code{diff}) between the coefficients 
   calculated from the synthetic and original data provides 
   a standardised \code{lack-of-fit = t(diff) \%*\% V^(-1) t(diff)}, where 
   \code{\%*\%} represents the matrix product and \code{V^(-1)} is the 
   inverse of the variance-covariance matrix for the mean coefficients 
   from the synthetic data. If the model used to synthesize the data is 
   correct this quantity, which is a Mahalanobis distance measure, will 
   follow a chi-squared distribution with degrees of freedom, and thus 
   expectation, equal to the number of parameters (\code{p}) in the fitted 
   model. When \code{incomplete.method = TRUE} the function above follows 
   a Hotelling's \code{T*2} distribution and the lack-of-fit statistic is 
   referred to an \code{F(p, m - p)}.}
 \item{lof.pvalue}{p-value for the combined lack-of-fit test of the NULL 
   hypothesis that the method used for synthesis retains all relationships 
   between variables that influence the parameters of the fit.}  
 \item{ci.plot}{\code{ggplot} of the the coefficients with confidence 
   intervals for models based on observed and synthetic data. 
   If \code{return.plot} was set to \code{FALSE} then \code{ci.plot} 
   is \code{NULL}.}
 \item{print.coef}{a logical value determining whether tables of 
   estimates for the original and synthetic data should be printed.}  
 \item{m}{the number of synthetic versions of the original (observed) data.} 
 \item{ncoef}{the number of coefficients in the fitted model (including an 
   intercept).} 
 \item{incomplete}{whether methods for incomplete synthesis due to Reiter 
   (2003) have been used in calculations.}
 \item{population.inference}{whether intervals as decribed by Karr et al. 
   (2016) have been calculated.}
} 

\references{
  Karr, A., Kohnen, C.N., Oganian, A., Reiter, J.P. and Sanil, A.P. (2006). 
  A framework for evaluating the utility of data altered to protect 
  confidentiality. \emph{The American Statistician}, \bold{60}(3), 224-232.

  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.
  
  Reiter, J.P. (2003) Inference for partially synthetic, public use microdata 
  sets. \emph{Survey Methodology}, \bold{29}, 181-188.
}

\seealso{
  \code{\link{summary.fit.synds}}
}

\examples{
ods <- SD2011[,c("sex","age","edu","smoke")]
s1 <- syn(ods, m = 3)
f1 <- glm.synds(smoke ~ sex + age + edu, data = s1, family = "binomial")
compare(f1, ods) 
compare(f1, ods, print.coef = TRUE, plot = "coef")
}
