\name{utility.gen}
\alias{utility.gen}
\alias{print.utility.gen}
\title{Distributional comparison of synthesised and observed data}
\description{
 Distributional comparison of synthesised data set with the original (observed) 
 data set using propensity scores.
}
\usage{
utility.gen(object, data, method = "logit", tree.method = "rpart", maxorder = 1,
            vars = NULL, aggregate = FALSE, maxit = 200, groups = FALSE, ngroups = 5,
            mincriterion = 0, nperms = 50, cp = 1e-8, minbucket = 5, digits = 2,
            print.zscores = FALSE, usethresh = TRUE, zthresh = 1.6, 
            print.variable.importance = FALSE, \dots)

\method{print}{utility.gen}(x, print.zscores = x$print.zscores, digits = x$digits, usethresh = x$usethresh, 
                            zthresh = x$zthresh, print.variable.importance = x$print.variable.importance, \dots)
}
\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s) as \code{object$syn}. This a single 
    data set when \code{object$m = 1} or a list of length \code{object$m} when 
    \code{object$m > 1}.}
  \item{data}{the original (observed) data set.}
  \item{method}{a single string or a vector of strings specifying the 
    method(s) for modeling the propensity scores.  Methods can be selected 
    from \code{"logit"} and \code{"cart"}. } 
  \item{tree.method}{Implementation of cart in R, can be \code{"rpart"} or 
    \code{"ctree"}.}
  \item{maxorder}{maximum order of interactions to be considered in 
    \code{"logit"} method. For model without interactions \code{0} should be
    provided.}  
  \item{vars}{variables to be included in the utility comparison. If none are 
    specified all the variables in the synthesised data will be included.}
  \item{aggregate}{logical flag as to whether the data should be aggregated by 
    collapsing identical rows before computation. This can lead to much faster 
    computation when all the variables are categorical. Only works for 
    \code{method = "logit"}.}
  \item{maxit}{maximum iterations to use when \code{method = "logit"}. If the 
    model does not converge in this number a warning will suggest increasing it.} 
  \item{groups}{logical variable as to whether any numerical variables are to 
    be categorised into groups.}
  \item{ngroups}{target number of groups for each variable: final number may 
    differ if there are many repeated values.}
  \item{mincriterion}{Criterion between 0 and 1 to use to control 
    \code{tree.method = "ctree"} when the tree will not be allowed to split 
    further. A value of \code{0.95} would be equivalent to a \code{5\%} significance test. 
    Here we set it to \code{0} to effectively disable this test and grow large trees.}
  \item{nperms}{number of permutations for the permutation test to obtain the 
    null distribution of the utility measure when \code{method = "cart"}.}
  \item{cp}{complexity parameter for classification with tree.method 
   \code{"rpart"}. Small values grow bigger trees.}
  \item{minbucket}{minimum number of observations allowed in a leaf for 
    classification when \code{method = "cart"}.}
  \item{\dots}{additional parameters passed to \code{\link{glm}},
    \code{\link{rpart}}, or \code{\link{ctree}}.} 
  \item{x}{an object of class \code{utility.gen}}
  \item{digits}{Number of digits to print in output}
  \item{print.zscores}{logical value as to whether z scores for coefficients of 
    the logit model should be printed.}
  \item{usethresh}{logical value as to whether printing of z scores below the 
    threshold should be suppressed.}
  \item{zthresh}{threshold value to use to suppress the printing of z scores 
    under \code{+}/\code{-} this value when \code{usethresh = TRUE} for 
    \code{method = "logit"}.}
  \item{print.variable.importance}{logical value as to whether the variable 
    importance measure should be printed when \code{tree.method = "rpart"}.}
}
\details{
  This function follows the method for evaluating the utility of masked data 
  as given in Snoke et al. (forthcoming) and originally proposed by Woo et al. 
  (2009). The original and synthetic data are combined into one dataset and 
  propensity scores, as detailed in Rosenbaum and Rubin (1983), are calculated 
  to estimate the probability of membership in the synthetic data set.
  The utility measure is based on the mean squared difference between these 
  probabilities and the probability expected if the data did not distinguish 
  the synthetic data from the original. The expected probability is just the 
  proportion of synthetic data in the combined data set, \code{0.5} when the 
  original and synthetic data have the same number of records.

  Propensity scores can be modeled by logistic regression \code{method = "logit"} 
  or by two different implementations of classification and regression trees as 
  method \code{"cart"}.  For logistic regression the predictors are all 
  variables in the data and their interactions up to order \code{maxorder}. The
  default of \code{1} gives all main effects and first order interactions. For logistic 
  regression the null distribution of the propensity score is derived and is 
  used to calculate ratios and standardised values.
  
  For \code{method = "cart"} the expectation and variance of the null 
  distribution is calculated from a permutation test.

  If missing values exist, indicator variables are added and included in the 
  model as recommended by Rosenbaum and Rubin (1984). For categorical variables, 
  \code{NA} is treated as a new category.}
  
\value{
  An object of class \code{utility.gen} which is a list including the utility 
  measures their expected null values for each synthetic set with the following 
  components:
  \item{call}{the call that produced the result.}
  \item{m}{number of synthetic data sets in object.}
  \item{method}{method used to fit propensity score.}
  \item{tree.method}{cart function used to fit propensity score when 
    \code{method = "cart"}.}
  \item{utilVal}{utility value(s).}
  \item{utilExp}{expected value(s) of the utility score if the synthesis method 
    is correct.}
  \item{utilR}{ratio(s) of \code{utilVal(s)} to \code{utilExp}.}
  \item{utilStd}{utility value standardised by expressing it as z scores, 
    difference(s) from the expected value divided by the expected standard 
    deviation.}
  \item{fit}{the fitted model for the propensity score or a list of fitted 
    models of length \code{m} if \code{m>0}.}
}
\references{
  Woo, M-J., Reiter, J.P., Oganian, A. and Karr, A.F. (2009). Global measures 
  of data utility for microdata masked for disclosure limitation. 
  \emph{Journal of Privacy and Confidentiality}, \bold{1}(1), 111-124.
  
  Rosenbaum, P.R. and Rubin, D.B. (1984). Reducing bias in observational studies 
  using subclassification on the propensity score.
  \emph{Journal of the American Statistical Association},
  \bold{79}(387), 516-524.
  
  Snoke, J., Raab, G.M., Nowok, B., Dibben, C. and Slavkovic, A. (forthcoming). 
  General and specific utility measures for synthetic data. 
}
\seealso{
 \code{\link{utility.tab}}
}
\examples{
  ods <- SD2011[1:1000, c("age", "bmi", "depress", "alcabuse", "englang")]
  s1 <- syn(ods, m = 5)
  utility.gen(s1, ods)
  u1 <- utility.gen(s1, ods)
  print(u1, print.zscores = TRUE, usethresh = TRUE)
  u2 <- utility.gen(s1, ods, groups = TRUE)
  print(u2, print.zscores = TRUE)
  u3 <- utility.gen(s1, ods, method = "cart", nperms = 20)
  print(u3, print.variable.importance = TRUE)
}
