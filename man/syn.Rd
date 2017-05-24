\name{syn}
\alias{syn}
\alias{syn.strata}
\alias{print.synds}
\title{Generating synthetic data sets}
\description{
  Generates synthetic version(s) of a data set. \code{syn.strata} performs 
  stratified synthesis.
}

\usage{
syn(data, method = vector("character", length = ncol(data)), 
    visit.sequence = (1:ncol(data)), predictor.matrix = NULL,  
    m = 1, k = nrow(data), proper = FALSE, minnumlevels = -1,  
    maxfaclevels = 60, rules = NULL, rvalues = NULL, 
    cont.na = NULL, semicont = NULL, smoothing = NULL, 
    event = NULL, denom = NULL, drop.not.used = FALSE, drop.pred.only = FALSE,
    default.method = c("normrank", "logreg", "polyreg", "polr"), models = FALSE,
    print.flag = TRUE, seed = "sample", ...)
    
syn.strata(data, strata = NULL, 
    minstratumsize = 10 + 10 * length(visit.sequence),  
    tab.strataobs = TRUE, tab.stratasyn = FALSE,
    method = vector("character", length = ncol(data)),
    visit.sequence = (1:ncol(data)), predictor.matrix = NULL,
    m = 1, k = nrow(data), proper = FALSE, minnumlevels = 5, 
    maxfaclevels = 60, rules = NULL, rvalues = NULL,
    cont.na = NULL, semicont = NULL, smoothing = NULL, 
    event = NULL, denom = NULL, drop.not.used = FALSE, drop.pred.only = FALSE,
    default.method = c("normrank","logreg","polyreg","polr"), models = FALSE,
    print.flag = TRUE, seed = "sample", ...)    
   
\method{print}{synds}(x, \dots)   
}

\arguments{
  \item{data}{a data frame or a matrix (\code{n} x \code{p}) containing the 
     original data. Observations are in rows and variables are in columns.}
  \item{method}{a single string or a vector of strings of length
     \code{ncol(data)} specifying the synthesising method to be
     used for each variable in the data. Order of variables is exactly the
     same as in \code{data}. If specified as a single string, the same method
     is used for all variables in a visit sequence unless a data type or
     a position in a visit sequence requires a different method. 
     If \code{method} is set to \code{"parametric"} the
     default synthesising method specified by the \code{default.method} argument
     are applied. Variables that are transformations of other variables can
     be synthesised using a passive method that is specified as a string 
     starting with \code{~} (see \code{\link{syn.passive}}). Variables that need 
     not to be synthesised have the empty method \code{""}. By default all 
     variables are synthesised using \code{"cart"} method, which is 
     \code{rpart} implementation of a CART model (see \code{\link{syn.cart}}). 
     See details for more information on method.}
  \item{visit.sequence}{a character vector of names of variables or an integer 
     vector of their column indices specifying the order of synthesis.
     The default sequence \code{1:ncol(data)} implies that column variables are
     synthesised from left to right. See details for more information.}
  \item{predictor.matrix}{a square matrix of size \code{ncol(data)} specifying
     the set of column predictors to be used for each target variable in the row.
     Each entry has value 0 or 1. A value of 1 means that the column
     variable is used as a predictor for the row variable. Order of
     variables is exactly the same as in \code{data}. By default all
     variables that are earlier in the visit sequence are used as predictors.
     For the default visit sequence (\code{1:ncol(data)}) the default
     \code{predictor.matrix} will have values of 1 in the lower triangle.
     See details for more information.}
  \item{m}{number of synthetic copies of the original (observed) data to be 
     generated. The default is \code{m = 1}.}
  \item{k}{a size of the synthetic data set (\code{k x p}),
     which can be smaller or greater than the size of the original data 
     set (\code{n x p}). The default is \code{nrow(data)} which means 
     that the number of individuals in the synthesised data is the same 
     as in the original (observed) data (\code{k = n}).}
  \item{proper}{a logical value with default set to \code{FALSE}. 
     If \code{TRUE} proper synthesis is conducted.}
  \item{minnumlevels}{a minimum number of values a numeric variable should
     have to be treated as numeric. Numeric variables with fewer levels
     than \code{minnumlevels} are changed into factors. If set to \code{-1}
     (default) numeric variables are left unchanged regardless of the number
     of values.}
  \item{maxfaclevels}{a maximum number of factor levels that can be handled.
     It can be increased but it may cause computational problems, especially 
     for parametric methods.}
  \item{rules}{a named list of rules for restricted values. Restricted values 
     are those that are determined explicitly by values of other variables.
     The names of the list elements must correspond to the variables names for 
     which the rules need to be specified.}
  \item{rvalues}{a named list of the values corresponding to the rules
     specified by \code{rules}.}
  \item{cont.na}{a named list of codes for missing values for continuous
     variables if different from the \code{R} missing data code \code{NA}.
     The names of the list elements must correspond to the variables names for 
     which the missing data codes need to be specified.}
  \item{semicont}{a named list of values at which semi-continuous
     variables have spikes. The names of the list elements must correspond 
     to the names of the semi-continuous variables.}     
  \item{smoothing}{a named list specifying smoothing method (\code{"density"} 
     or \code{""}) to be used for selected variables. Smoothing can only be 
     applied to continuous variables synthesised using \code{sample}, 
     \code{ctree}, \code{cart} or \code{normrank} method. The names of the 
     list elements must correspond to the names of the variables whose values 
     are to be smoothed. Smoothing is applied to the synthesised values.
     For \code{"density"} smoothing a Gaussian kernel density estimator is 
     applied with bandwidth selected using the Sheather-Jones 
     'solve-the-equation' method (see \code{\link{bw.SJ}}).}  
  \item{event}{a named list specifying for survival data the names of 
     corresponding event indicators. The names of the list elements 
     must correspond to the names of the survival variables.}
  \item{denom}{a named list specifying for variables to be modelled using 
     binomial regression the names of corresponding denominator variables. 
     The names of the list elements must correspond to the names of the 
     variables to be modelled using binomial regression.}
  \item{drop.not.used}{a logical value. If \code{TRUE} (default) variables not 
     used in synthesis are not saved in the synthesised data and are not 
     included in the corresponding synthesis parameters.}
  \item{drop.pred.only}{a logical value. If \code{TRUE} (default) variables not 
     synthesised and used as predictors only are not saved in the synthesised 
     data.}
  \item{default.method}{a vector of four strings containing the default
     parametric synthesising methods for numerical variables, factors 
     with two levels, unordered factors with more than two levels 
     and ordered factors with more than two levels respectively. 
     They are used when \code{method} is set to \code{"parametric"} or 
     when there is an inconsistency between variable type and provided method.}
  \item{models}{if \code{TRUE} estimates of models fitted to the original data 
     and used to generate the synthetic values are stored.}
  \item{print.flag}{if \code{TRUE} (default) synthesising history and 
     information messages will be printed at the console. For silent 
     computation use \code{print.flag = FALSE}.}
  \item{seed}{an integer to be used as an argument for the \code{set.seed()}.
     If no integer is provided, the default \code{"sample"} will generate one
     and it will be stored. To prevent generating an integer set \code{seed}
     to \code{NA}.}
  \item{\dots}{additional arguments to be passed to synthesising functions. See 
     section 'Details' below for more information.}
  \item{strata}{a numeric vector with strata identifiers or a string vector with 
     names of stratifying variable(s).} 
  \item{minstratumsize}{minimum size of each stratum.}    
  \item{tab.strataobs}{a logical value indicating whether a frequency table of
     the number of observations in strata in the original data set should be 
     printed.}   
  \item{tab.stratasyn}{a logical value indicating whether a frequency table of
     the number of observations in strata in the synthetic data set(s) should 
     be printed.}           
  \item{x}{an object of class \code{synds}; a result of a call to \code{syn}.}
}

\details{
 Only variables that are in \code{visit.sequence} with corresponding non-empty 
 \code{method} are synthesised. The only exceptions are event indicators. They 
 are synthesised along with the corresponding time to event variables and should
 not be included in \code{visit.sequence}. All other variables (not in 
 \code{visit.sequence} or in \code{visit.sequence} with a corresponding blank 
 method) can be used as predictors. Including them in \code{visit.sequence} 
 generates a default \code{predictor.matrix} reflecting the order of variables 
 in the \code{visit.sequence} otherwise \code{predictor.matrix} has to be 
 adjusted accordingly. All predictors of the variables that are not in
 \code{visit.sequence} or are in \code{visit.sequence} but with a blank method
 are removed from \code{predictor.matrix}.

 Variables to be synthesised that are not synthesised yet cannot be used
 as predictors. Also all variables used in passive synthesis or in restricted
 values rules (\code{rules}) have to be synthesised before the variables they
 apply to.

 Mismatch between data type and synthesising method stops execution and
 print an error message but numeric variables with number of levels less
 than \code{minnumlevels} are changed into factors and methods are changed
 automatically, if necessary, to methods for categorical variables. 
 Methods for variables not in a visit sequence will be changed into blank.

 The built-in elementary synthesising methods include:
 \describe{
   \item{ctree, cart}{classification and regression trees (CART)}
   \item{survctree}{classification and regression trees (CART)
      for duration time data (parametric methods for survival data are 
      not implemented yet)}
   \item{norm}{normal linear regression}
   \item{normrank}{normal linear regression preserving the marginal distribution}
   \item{lognorm, sqrtnorm, cubertnorm}{normal linear regression after 
      natural logarithmic, square root and cube root transformation of 
      a dependent variable respectively} 
   \item{logreg}{logistic regression}
   \item{polyreg}{unordered polytomous regression}
   \item{polr}{ordered polytomous regression}
   \item{pmm}{predictive mean matching}
   \item{sample}{random sample from the observed data}
   \item{passive}{function of other synthesised data}
   \item{nested}{bootstrap sample within each category of the original grouping 
     variable)}
  }
 The functions corresponding to these methods are called \code{syn.method}, 
 where \code{method} is a string with the name of a synthesising method. 
 For instance a function corresponding to \code{ctree} function is called
 \code{syn.ctree}. A new synthesising method can be introduced by writing 
 a function named \code{syn.newmethod} and then specifying \code{method} 
 parameter of \code{syn} function as \code{"newmethod"}.
 
 In order to use \code{"nested"} sampling, \code{method} parameter of \code{syn} 
 function has to be specified as \code{"nested.varname"}, where \code{"varname"}
 is the name of the grouped (less detailed) variable, the only one used in 
 nested synthesis. A variable synthesised using \code{"nested"} method is 
 excluded from synthesising other variables except when used for \code{"nested"} 
 method.  
 
 Additional parameters can be passed to synthesising methods as part of the 
 \code{dots} argument. They have to be named using period-separated method and 
 parameter name (\code{method.parameter}). For instance, in order to set 
 a \code{minbucket} (minimum number of observations in any terminal node of 
 a CART model) for a \code{ctree} synthesising method, \code{ctree.minbucket}
 has to be specified. The parameters are method-specific and will be used for 
 all variables to be synthesised using that method. See help for 
 \code{syn.method} for further details about the allowed parameters for 
 a specific method.  
}

\value{
  An object of class \code{synds}, which stands for 'synthesised 
  data set'. It is a list with the following components:

  \item{call}{an original call to \code{syn}.}
  \item{m}{number of synthetic versions of the original (observed) data.}
  \item{syn}{a data frame (for \code{m = 1}) or a list of \code{m} data frames 
    (for \code{m > 1}) with synthetic data set(s).}
  \item{method}{a vector of synthesising methods applied to each variable
     in the saved synthesised data.}
  \item{visit.sequence}{a vector of column indices of the visiting sequence.
     The indices refer to the columns in the saved synthesised data.}
  \item{predictor.matrix}{a matrix specifying the set of predictors used
     for each variable in the saved synthesised data.}
  \item{smoothing}{a vector specifying smoothing methods 
     applied to each variable in the saved synthesised data.}
  \item{event}{a vector of integers specifying for survival data the column 
     indices for corresponding event indicators. The indices refer to the 
     columns in the saved synthesised data.}
  \item{denom}{a vector of integers specifying for variables modelled 
     using binomial regression the column indices for corresponding 
     denominator variables. The indices refer to the columns in the 
     saved synthesised data.}
  \item{proper}{a logical value indicating whether proper synthesis was
     conducted.}
  \item{n}{a number of cases in the original data.}  
  \item{k}{a number of cases in the synthesised data.}
  \item{rules}{a list of rules for restricted values applied to the synthetic 
    data.}
  \item{rvalues}{a list of the values corresponding to the rules
     specified by \code{rules}.}
  \item{cont.na}{a list of codes for missing values for continuous variables.}
  \item{semicont}{a list of values for semi-continuous variables at which they 
     have spikes.}
  \item{drop.not.used}{a logical value indicating whether variables not used
     in synthesis are saved in the synthesised data and corresponding
     synthesis parameters.}
  \item{drop.pred.only}{a logical value indicating whether variables not
     synthesised and used as predictors only are saved in the synthesised data.}
  \item{seed}{an integer used as a \code{set.seed()} argument.}   
  \item{var.lab}{a vector of variable labels for data imported from SPSS using 
     \code{read.obs()}.}   
  \item{val.lab}{a list value labels for factors for data imported from SPSS 
     using \code{read.obs()}.}       
  \item{obs.vars}{a vector of all variable names in the observed data set.}
  \item{models}{estimates of models fitted to the original data and used to 
     generate the synthetic values.}
}
\note{See package vignette for additional information.}

\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.
}

\seealso{\code{\link{compare.synds}}, \code{\link{summary.synds}}}

\examples{
### selection of variables
vars <- c("sex","age","marital","income","ls","smoke")
ods  <- SD2011[1:2000,vars]
 
### default synthesis
s1 <- syn(ods)
s1
  
### synthesis with default parametric methods
s2 <- syn(ods, method = "parametric", seed = 1)
s2$method
  
### multiple synthesis of selected variables with customised methods
s3 <- syn(ods, visit.sequence = c(2, 1, 4, 5), m = 2,
          method = c("logreg","sample","","normrank","ctree",""),
          ctree.minbucket = 10)
summary(s3)
summary(s3, msel = 1:2)
  
### adjustment to the default predictor matrix 
s4.ini <- syn(data = ods, visit.sequence = c(1, 2, 5, 3),
              m = 0, drop.not.used = FALSE)
pM.cor <- s4.ini$predictor.matrix
pM.cor["marital","ls"] <- 0
s4 <- syn(data = ods, visit.sequence = c(1, 2, 5, 3),
          predictor.matrix = pM.cor)
  
### handling missing values in continuous variables
s5 <- syn(ods, cont.na = list(income = c(NA, -8)))
  
### rules for restricted values - marital status of males under 18 should be 'single'
s6 <- syn(ods, rules = list(marital = "age < 18 & sex == 'MALE'"),
          rvalues = list(marital = 'SINGLE'), method = "parametric", seed = 1)
with(s6$syn, table(marital[age < 18 & sex == 'MALE']))
### results for default parametric synthesis without the rule  
with(s2$syn, table(marital[age < 18 & sex == 'MALE']))

### stratified synthesis
s7 <- syn.strata(ods, strata = "sex")
}

\keyword{datagen}
\keyword{tree}
\keyword{regression}
