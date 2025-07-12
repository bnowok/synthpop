\name{utility.gen}
\alias{utility.gen}
\alias{utility.gen.synds}
\alias{utility.gen.data.frame}
\alias{utility.gen.list}
\alias{print.utility.gen}
\title{Distributional comparison of synthesised and observed data}
\description{
 Distributional comparison of synthesised data set with the original (observed) 
 data set using propensity scores.
 
 This function can be also used with synthetic data NOT created by 
 \code{syn()}, but then additional parameters \code{not.synthesised} 
 and \code{cont.na} might need to be provided.
}
\usage{
\method{utility.gen}{synds}(object, data, 
            method = "cart", maxorder = 1, k.syn = FALSE, tree.method = "rpart",
            max.params = 400, print.stats = c("pMSE", "S_pMSE"), resamp.method = NULL, 
            nperms = 50, cp = 1e-3, minbucket = 5, mincriterion = 0, vars = NULL, 
            aggregate = FALSE, maxit = 200, ngroups = NULL, print.flag = TRUE,
            print.every = 10, digits = 6, print.zscores = FALSE, zthresh = 1.6,
            print.ind.results = FALSE, print.variable.importance = FALSE, \dots)

\method{utility.gen}{data.frame}(object, data, not.synthesised = NULL, cont.na = NULL, 
            method = "cart", maxorder = 1, k.syn = FALSE, tree.method = "rpart",
            max.params = 400, print.stats = c("pMSE", "S_pMSE"), resamp.method = NULL, 
            nperms = 50, cp = 1e-3, minbucket = 5, mincriterion = 0, vars = NULL, 
            aggregate = FALSE, maxit = 200, ngroups = NULL, print.flag = TRUE,
            print.every = 10, digits = 6, print.zscores = FALSE, zthresh = 1.6,
            print.ind.results = FALSE, print.variable.importance = FALSE, \dots)

\method{utility.gen}{list}(object, data, not.synthesised = NULL, cont.na = NULL, 
            method = "cart", maxorder = 1, k.syn = FALSE, tree.method = "rpart",
            max.params = 400, print.stats = c("pMSE", "S_pMSE"), resamp.method = NULL, 
            nperms = 50, cp = 1e-3, minbucket = 5, mincriterion = 0, vars = NULL, 
            aggregate = FALSE, maxit = 200, ngroups = NULL, print.flag = TRUE,
            print.every = 10, digits = 6, print.zscores = FALSE, zthresh = 1.6,
            print.ind.results = FALSE, print.variable.importance = FALSE, \dots)


\method{print}{utility.gen}(x, digits = NULL, zthresh = NULL, 
               print.zscores = NULL, print.stats = NULL,
               print.ind.results = NULL, print.variable.importance = NULL, \dots)
}

\arguments{
  \item{object}{it can be an object of class \code{synds}, which stands for 
    'synthesised data set'. It is typically created by function \code{syn()} 
    and it includes \code{object$m} synthesised data set(s) as \code{object$syn}. 
    This a single data set when \code{object$m = 1} or a list of length 
    \code{object$m} when \code{object$m > 1}. Alternatively, when data are 
    synthesised not using \code{syn()}, it can be a data frame with a synthetic 
    data set or a list of data frames with synthetic data sets, all created from 
    the same original data with the same variables and the same method.}
  \item{data}{the original (observed) data set.}
  \item{not.synthesised}{a vector of variable names for any variables that has 
    been left unchanged in the synthetic data. Not required if oject is of 
    class \code{synds}}
  \item{cont.na}{a named list of codes for missing values for continuous
    variables if different from the \code{R} missing data code \code{NA}.
    The names of the list elements must correspond to the variables names for 
    which the missing data codes need to be specified. Not required if oject is
    of class \code{synds}}  
  \item{method}{a single string specifying the method for modeling the propensity 
    scores. Method can be selected from \code{"logit"} and \code{"cart"}. } 
  \item{maxorder}{maximum order of interactions to be considered in 
    \code{"logit"} method. For model without interactions \code{0} should be
    provided.}   
  \item{k.syn}{a logical indicator as to whether the sample size itself has 
    been synthesised.}
  \item{tree.method}{implementation of \code{"cart"} method that is used when 
    \code{method = "cart"}. It can be \code{"rpart"} or \code{"ctree"}.}
  \item{max.params}{the maximum number of parameters for a \code{"logit"} model 
    which alerts the user to possible fitting failure.}
  \item{print.stats}{statistics to be printed must be a selection from
    \code{"pMSE"}, \code{"SPECKS"}, \code{"PO50"}, \code{"S_pMSE"}, 
    \code{"S_SPECKS"}, \code{"S_PO50"}. If \code{print.stats = "all"}, 
    all of the measures mentioned above will be printed.}
  \item{resamp.method}{method used for resampling estimates of standardized 
    measures can be \code{"perm"}, \code{"pairs"} or \code{"none"}.  
    Defaults to \code{"pairs"} if \code{print.stats} includes \code{"S_SPECKS"} 
    or \code{"S_PO50"} or synthesis is incomplete else defaults to \code{"perm"}
    if method is \code{"cart"} or to \code{NULL}, no resampling needed, 
    if method is \code{"logit"}. \code{"none"} can be used to get results 
    without standardized measures e.g. in simulations.}
  \item{nperms}{number of permutations for the permutation test to obtain the 
    null distribution of the utility measure when \code{resamp.method = "perm"}.}
  \item{cp}{complexity parameter for classification with tree.method 
    \code{"rpart"}. Small values grow bigger trees.}
  \item{minbucket}{minimum number of observations allowed in a leaf for 
    classification when \code{method = "cart"}.}
  \item{mincriterion}{criterion between 0 and 1 to use to control 
    \code{tree.method = "ctree"} when the tree will not be allowed to split 
    further. A value of \code{0.95} would be equivalent to a \code{5\%} 
    significance test. Here we set it to \code{0} to effectively disable this 
    test and grow large trees.}
  \item{vars}{variables to be included in the utility comparison. It can be 
    a character vector of names of variables or an integer vector of their 
    column indices. If none are specified all the variables in the synthesised 
    data will be included.}
  \item{aggregate}{logical flag as to whether the data should be aggregated by 
    collapsing identical rows before computation. This can lead to much faster 
    computation when all the variables are categorical. Only works for 
    \code{method = "logit"}.}
  \item{maxit}{maximum iterations to use when \code{method = "logit"}. If the 
    model does not converge in this number a warning will suggest increasing 
    it.} 
  \item{ngroups}{target number of groups for categorisation of each numeric 
    variable: final number may differ if there are many repeated values. If 
    \code{NULL} (default) variables are not categorised into groups.}
  \item{print.flag}{TRUE/FALSE to indicate if any messages should 
    be printed during calculations. Change to FALSE for simulations.}  
  \item{print.every}{controls the printing of progress of resampling when  
    \code{resamp.method} is not \code{NULL}. When \code{print.every = 0} 
    no progress is reported, otherwise the resample number is printed every 
    \code{print.every}.}
  \item{\dots}{additional parameters passed to \code{\link[stats]{glm}},
    \code{\link[rpart]{rpart}}, or \code{\link[party]{ctree}}.} 
  \item{x}{an object of class \code{utility.gen}.}
  \item{digits}{number of digits to print in the default output values.}  
  \item{zthresh}{threshold value to use to suppress the printing of z-scores 
    under \code{+}/\code{-} this value for \code{method = "logit"}. If set to 
    \code{NA} all z-scores are printed.}
  \item{print.zscores}{logical value as to whether z-scores for coefficients of 
    the logit model should be printed.}  
  \item{print.ind.results}{logical value as to whether utility score results 
    from individual syntheses should be printed.}
  \item{print.variable.importance}{logical value as to whether the variable 
    importance measure should be printed when \code{tree.method = "rpart"}.}
}
\details{
  This function follows the method for evaluating the utility of masked data 
  as given in Snoke et al. (2018) and originally proposed by Woo et al. 
  (2009). The original and synthetic data are combined into one dataset and 
  propensity scores, as detailed in Rosenbaum and Rubin (1983), are calculated 
  to estimate the probability of membership in the synthetic data set.
  The utility measure is based on the mean squared difference between these 
  probabilities and the probability expected if the data did not distinguish 
  the synthetic data from the original. 
 
  If \code{k.syn = FALSE} the expected probability is just the proportion of 
  synthetic data in the combined data set, \code{0.5} when the original and 
  synthetic data have the same number of records. Setting \code{k.syn = TRUE} 
  indicates that the numbers of observations in the synthetic data was 
  synthesised and not fixed by the synthesiser. In this case the expected 
  probability will be \code{0.5} in all cases and the model to discriminate
  between observed and synthetic will include an intercept term. This will 
  usually only apply when the standalone version of this function 
  \code{utility.gen.sa()} is used.

  Propensity scores can be modeled by logistic regression \code{method = "logit"} 
  or by two different implementations of classification and regression trees as 
  method \code{"cart"}. For logistic regression the predictors are all variables 
  in the data and their interactions up to order \code{maxorder}. The default of 
  \code{1} gives all main effects and first order interactions. For logistic 
  regression the null distribution of the propensity score is derived and is 
  used to calculate ratios and standardised values.
  
  For \code{method = "cart"} the expectation and variance of the null 
  distribution is calculated from a permutation test. Our recent work
  indicates that this method can sometimes give misleading results.

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
  \item{resamp.method}{type of resampling used to get \code{pMSEExp} and 
    \code{pval}.} 
  \item{maxorder}{see above.}
  \item{vars}{see above.}
  \item{nfix}{see above.}
  \item{aggregate}{see above.}
  \item{maxit}{see above.}
  \item{ngroups}{see above.}
  \item{df}{degrees of freedom for the chi-squared test for logit models
    derived from the number of non-aliased coefficients in the logistic model, 
    minus \code{1} for \code{k.syn = FALSE}.}
  \item{mincriterion }{see above.}
  \item{nperms}{see above.}
  \item{incomplete}{TRUE/FALSE indicator if any of the variables being compared 
    are not synthesised.}
  \item{pMSE}{propensity score mean square error from the utility model or a  
    vector of these values if \code{object$m > 1}.}
  \item{S_pMSE}{ratio(s) of \code{pMSE} to its Null expectation.}
  \item{PO50}{percentage  over 50\% of each synthetic data set where the model 
    used correctly predicts whether real or synthetic.}
  \item{S_PO50}{ratio(s) of \code{PO50} to its Null expectation.}
  \item{SPECKS}{Kolmogorov-Smirnov statistic to compare the propensity scores  
    for the original and synthetic records.}
  \item{S_SPECKS}{ratio(s) of \code{SPECKS} to its Null expectation.}
  \item{print.stats}{see above.}  
  \item{fit}{the fitted model for the propensity score or a list of fitted 
    models of length \code{m} if \code{m > 0}.}
  \item{nosplits}{for resampling methods and cart models, a list of the number
    of times from the total each resampled cart model failed to select any 
    splits to classify the indicator. Indicates that this method is not working 
    correctly and results should not be used but a logit model selected instead.}
  \item{digits}{see above.}
  \item{print.ind.results}{see above.}
  \item{print.zscores}{see above.}
  \item{zthresh}{see above.}
  \item{print.variable.importance}{see above.}
}
\references{
  Woo, M-J., Reiter, J.P., Oganian, A. and Karr, A.F. (2009). Global measures 
  of data utility for microdata masked for disclosure limitation. 
  \emph{Journal of Privacy and Confidentiality}, \bold{1}(1), 111-124.
  
  Rosenbaum, P.R. and Rubin, D.B. (1984). Reducing bias in observational 
  studies using subclassification on the propensity score.
  \emph{Journal of the American Statistical Association},
  \bold{79}(387), 516-524.
  
  Snoke, J., Raab, G.M., Nowok, B., Dibben, C. and Slavkovic, A. (2018). 
  General and specific utility measures for synthetic data. 
  \emph{Journal of the Royal Statistical Society: Series A}, \bold{181}, 
  Part 3, 663-688.
}
\seealso{
  \code{\link{utility.tab}}    
}
\examples{
\dontrun{
  ods <- SD2011[1:1000, c("age", "bmi", "depress", "alcabuse", "nofriend")]
  s1 <- syn(ods, m = 5, method = "parametric", 
            cont.na = list(nofriend = -8))
    
  ### synthetic data provided as a 'synds' object   
  u1 <- utility.gen(s1, ods)
  print(u1, print.zscores = TRUE, zthresh = 1, digits = 6)
  u2 <- utility.gen(s1, ods, ngroups = 3, print.flag = FALSE)
  print(u2, print.zscores = TRUE)
  u3 <- utility.gen(s1, ods, method = "cart", nperms = 20)
  print(u3, print.variable.importance = TRUE)
    
  ### synthetic data provided as 'list'
  utility.gen(s1$syn, ods, cont.na = list(nofriend = -8))  
  }
}
