\name{disclosure}
\alias{disclosure}
\alias{disclosure.synds}
\alias{disclosure.data.frame}
\alias{disclosure.list}
\alias{print.disclosure}
\title{Disclosure measures}
\description{
  Calculates disclosure measures for synthetic data. 
  NOTE: The other function that calculates disclosure results for multiple
  targets has been renamed as multi.disclosure from disclosure.summary.
}
\usage{
\method{disclosure}{synds}(object, data, keys , target , print.flag = TRUE,
           denom_lim = 5, exclude_ov_denom_lim = FALSE, not.targetlev = NULL,
           usetargetNA = TRUE, usekeysNA = TRUE, 
           exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
           ngroups_target = NULL, ngroups_keys = NULL, 
           thresh_1way = c(50, 90),thresh_2way = c(4, 80),
           digits = 2, to.print =c("short"),\dots) 

\method{disclosure}{data.frame}(object, data,cont.na = NULL, keys , target , print.flag = TRUE,
           denom_lim = 5, exclude_ov_denom_lim = FALSE, 
           not.targetlev = NULL,
           usetargetNA = TRUE, usekeysNA = TRUE, 
           exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
           ngroups_target = NULL, ngroups_keys = NULL, 
           thresh_1way = c(50, 90),thresh_2way = c(4, 80),
           digits = 2, to.print =c("short"), compare.synorig = TRUE, \dots) 

\method{disclosure}{list}(object, data,cont.na = NULL, keys , target , print.flag = TRUE,
           denom_lim = 5, exclude_ov_denom_lim = FALSE, 
           not.targetlev = NULL,
           usetargetNA = TRUE, usekeysNA = TRUE, 
           exclude.keys =NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
           ngroups_target = NULL, ngroups_keys = NULL, 
           thresh_1way = c(50, 90),thresh_2way = c(4, 80),
           digits = 2, to.print =c("short"), compare.synorig = TRUE, \dots) 
           
\method{print}{disclosure}(x,  to.print =NULL, digits = NULL, \dots)
}

\arguments{
  \item{object}{an object of class \code{synds}, which stands for 
  'synthesised data set'. It is typically created by function \code{syn()} and 
  it includes  \code{object$m} synthesised data set(s) as \code{object$syn}. 
  This a single data set when \code{object$m = 1} or a list of length 
  \code{object$m} when \code{object$m > 1}. Alternatively, when data are 
  synthesised not using \code{syn()}, it can be a data frame with a synthetic 
  data set or a list  of data frames with synthetic data sets, all created from 
  the same original data with the same variables.}  
  \item{data}{the original (observed) data set.}
  \item{cont.na}{For data NOT supplied as a synthetic data object created by 
  \code{synthpop}, this gives special values for continuous variables as
   described in the documentation for the function \code{syn}.}
  \item{keys}{vector of variable names or  column numbers in data 
  that are also present in the synthetic data to act as 
       quasi-identifiers for identity or attribute disclosure.}
  \item{target}{name of target variable for attribute disclosure.}
  \item{denom_lim}{Limit to use to exclude large key-target group, see next item.}
  \item{exclude_ov_denom_lim}{logical to exclude key targetcombinations
    that contribute more than  \code{denom_lim} disclosive records.
    These are often flagged from  \code{thresh_2way}  where the 
    first element corresponds to  \code{denom_lim}}
  \item{print.flag}{logical value as to whether a line is printed as
    disclosure is calculated for each synthetic data set.}
  \item{digits}{number of digits to print for disclosure measures.} 
  \item{usetargetNA}{determines whether NA values in target are to be used 
            in checking for disclosure}
  \item{usekeysNA}{determines whether NA values in keys are to be used 
            in checking for disclosure.}
  \item{not.targetlev}{Character variable giving level of target to be
   excluded from disclosure measures. Usually identified by checklev_1way.}
  \item{exclude.keys}{vector of names of keys that, with the next two items
     will define the target and key combinations to be excluded from the
     calculation of disclosure measures. Often identified by checklev_2way.}
  \item{exclude.keylevs}{vector of the same length as exclude.keys that give 
      the levels  to be excluded for the corresponding key.}
  \item{exclude.targetlevs}{vector of the same length as exclude.keys that give 
      the levels  of target that will be excluded for each key and key level.}
  \item{ngroups_target}{Unless set to NULL (the default) a numeric target variable
     will be grouped into \code{ngroups_target} categories.}
  \item{ngroups_keys}{Unless set to NULL (the default) any numeric variable
     will be grouped into categories.  If \code{ngroups_keys} is of length 1 all numeric
     keys will be have the same number of groups. Otherwise \code{ngroups_keys}
     needs to be the same length as keys and will give the number of groups for each
     key. If an element of \code{ngroups_keys} is zero, no grouping will be done.}
  \item{thresh_1way}{A vector of two numeric values both of which meed to be exceeded
   for warnings about a level of the target that may be dominating the results. 
   The first is the count of all disclosive records for this level of the target,
   and the second is the \% of all original records for this level of the target. Default is c(50, 90), 
   meaning a group of 50 disclosive records for this level of the target
   where they make up over 90\% of all disclosive records.}
  \item{thresh_2way}{A vector of two numeric values both of which meed to be 
   exceeded for warnings about a level of the target that may be dominating 
   the results. The first is the count of disclosive records for a 
   quasi-identifier used to identify possible s that
   are searched for the most disclosive key-target combination. The second is the  
   percentage of all original records for each combination examined that must
   be exceeded to trigger a warning.  Default is c(5, 80), 
   meaning a pairs found from key-target groups of more than 5 records 
   where over 80\% of all the original values with these key-target 
   pairs have this level of the target.}
\item{to.print}{Vector to determine what aspects of an object of
   class disclosure will be printed. Must consist  of  one or more of the 
   following "short", "ident", "attrib","allCAPs", "all", "check_1way", 
   "check_2way", "exclusions". Default is "short" giving a brief summary.}
  \item{compare.synorig}{a logical value  to determine if the functions
  \code{synorig.compare()} should be used to check that data sets can be 
   compared. Used when the synthetic data are supplied as a data.frame or 
   a list when  default set to TRUE.} 
\item{\dots}{additional parameters}  
 \item{x}{an object of class \code{disclosure}.}
}

\details{ Calculates identity disclosure measures for a for a set of keys,
         (quasi identifiers) and attribute disclosure measures for  one
         variable from the same  set of keys considered as a target. The 
         function \code{\link{multi.disclosure}} calls this function and 
         summarises the attribute disclosure measures for multiple targets.
         See the vignette 
}

\value{An object of class \code{disclosure} which is a list with the following 
  components. 
  \item{call}{the call that created the object.}
  \item{ident}{Table of measures of identity disclosure
              one for each synthesis. Measures are "UiO","UiS","UiSiO" and
              "repU". See vignette disclosure.pdf for an explanation of these 
              and the following measures.}
  \item{attrib}{Table of measures of attribute disclosure  one for each 
    synthesis. These include "DiO","DiS","iSO","DiSCO" and "DiSDiO".
    The measures "DiO" and "DiS" are the percentage of the target that are
    disclosed from the original and synthetic data with these keys. The next
    measure "iSO" gives the percentage of the  key combinations 
    in the synthetic data that are present in the original - 
    one was in which the disclosure.  "DiSCO" gives the percentage 
    of original records where the attribution to the target is correct as judged 
    from the original. "DiSDiO" gives the \% of origina; records in "DISCO"
    that are unique in the original data.
    The table also as gives the maximum and mean of the denominators for the
    "DiSCO" measure i.e. the distribution for every record that leads to a 
    correct disclosure of the number of observations with
    the same keys and the same correct target in the synthetic data. 
    Large denominators are often an   indication that the disclosure is 
    something that might be expected from   prior knowledge of relations.} 
  \item{allCAPs}{Table of the following measures of correct attribution
      probability: "baseCAPd","CAPd", "CAPs" , "DCAP" and  "TCAP"'} 
  \item{check_1way}{A data frame with one record per synthesis
      identifying the level of the target with numbers of disclosive records
      that are above thresholds defined by \code{thresh_1way}, with default
      value c(50,90). This means that there must be more than 50 disclosive       records
    with this level of the target, and that 90\% or more of all disclosive 
    records must have this target. The value of most_dis_lev will be blank 
    if no level exceeds these thresholds.
    Note this level will be identified  for data without 
    excluded or missing values of keys if there are any excluded records.} 
  \item{check1}{The level of the target identified by check_1way
   ` or blank if none}
  \item{check_2way}{A list of length number of syntheses giving details 
    for each of the two-way combinations of target and keys where the
    the numbers of disclosive records are above thresholds defined by 
    \code{thresh_2way}. The default value for this is c(5, 80),
    meaning that there must be at least 5 records with this combination
    of targets and keys and that 80\% or more of records in the original
    data with this level of the key will have this level of the target.
    If no combinations exceed \code{thresh_2way} for one of the syntheses 
    then the list element is NULL. 
    Such disclosive combinations are often associated with a high prior 
    probability of   the target from just this level of one of the  keys 
    in the original data. 
    Note these combinations will be identified  for data without 
    excluded or missing values of keys if there are any excluded combinations 
    or target if any of \code{usekeysNA} or \code{usetargetNA} are FALSE.}
  \item{Nexclusions}{A list of length number of syntheses with number of 
      records excluded from attribute measures for different reasons.}
  \item{keys}{as input}
  \item{digits}{as input}
  \item{Norig}{Number of records in data}
  \item{to.print}{as input}
  }
  
\note{See package vignette disclosure.pdf for additional information including
  formal definitions of all quantities and worked examples.}

\references{
 See references in package vignette}
 
\seealso{  \code{\link{syn}} \code{\link{multi.disclosure}}
}

\examples{
library(synthpop)
ods <- SD2011[, c("sex", "age", "edu", "marital", "income")]
odsF <- numtocat.syn(ods, numtocat = "income", catgroups = 7, cont.na = list(income = -8))
s1 <- syn(odsF$data, method = "ctree",seed = 75, m=3, k=1000)
disc1 <- disclosure(s1, odsF$data, target = "income", 
                    keys = c("sex", "age", "edu","marital"))
}
