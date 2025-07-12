\name{multi.disclosure}
\alias{multi.disclosure}
\alias{multi.disclosure.synds}
\alias{multi.disclosure.data.frame}
\alias{multi.disclosure.list}
\alias{print.multi.disclosure}
\title{Disclosure measures for multiple of target variables.}
\description{
  Calculates, prints and plots tables of disclosure measures for a set of
  target variables from a fixed set of keys to form quasi-identifiers. 
  The calculations of disclosure measures are done by the function
  \code{\link{disclosure}} for each target.

  This function can be also used with synthetic data NOT created by 
  \code{syn()}, or even made anonymous by other methods such as sampling
  More details of the measures calculated can be found in the package vignette
  "Disclosure measures for Synthetic Data".
}

\usage{
\method{multi.disclosure}{synds}(object, data, 
           keys , targets = NULL, print.flag = TRUE, 
           denom_lim = 5, exclude_ov_denom_lim = FALSE,
           not.targetslev = NULL,  
           usetargetsNA = TRUE,  usekeysNA = TRUE, 
           exclude.keys = NULL, exclude.keylevs = NULL,  exclude.targetlevs = NULL,
           ngroups_targets = NULL, ngroups_keys = NULL, 
           ident.meas = "repU", attrib.meas = "DiSCO",
           thresh_1way = c(50, 90),thresh_2way = c(4, 80), 
           digits = 2, plot = TRUE,  \dots)

               
\method{multi.disclosure}{data.frame}(object, data, cont.na = NULL, 
           keys , targets = NULL,  print.flag = TRUE, 
           denom_lim = 5, exclude_ov_denom_lim = FALSE,
           not.targetslev = NULL, 
           usetargetsNA = TRUE,  usekeysNA = TRUE, 
           exclude.keys = NULL, exclude.keylevs = NULL,  exclude.targetlevs = NULL,
           ngroups_targets = NULL, ngroups_keys = NULL, 
           ident.meas = "repU", attrib.meas = "DiSCO",
           thresh_1way = c(50, 90),thresh_2way = c(4, 80), 
           digits = 2, plot = TRUE,  compare.synorig = TRUE,  \dots)

\method{multi.disclosure}{list}(object, data, cont.na = NULL,
            keys , targets = NULL,  print.flag = TRUE, 
            denom_lim = 5, exclude_ov_denom_lim = FALSE,
           not.targetslev = NULL,  
           usetargetsNA = TRUE,  usekeysNA = TRUE, 
           exclude.keys = NULL, exclude.keylevs = NULL, exclude.targetlevs = NULL,
           ngroups_targets = NULL, ngroups_keys = NULL, 
           ident.meas = "repU", attrib.meas = "DiSCO",
           thresh_1way = c(50, 90),thresh_2way = c(4, 80), 
           digits = 2, plot = TRUE, compare.synorig = TRUE,\dots)


\method{print}{multi.disclosure}(x, digits = NULL, plot = NULL, to.print =  c("ident","attrib"),
       \dots)}


\arguments{
  \item{object}{an object of class \code{synds}, which stands for 'synthesised 
    data set'. It is typically created by function \code{syn()} and it includes 
    \code{object$m} synthesised data set(s) as \code{object$syn}. This a single 
    data set when \code{object$m = 1} or a list of length \code{object$m} when 
    \code{object$m > 1}. Alternatively, when data are synthesised not using 
    \code{syn()}, it can be a data frame with a synthetic data set or a list 
    of data frames with synthetic data sets, all created from the same original 
    data with the same variables and the same method.}  
  \item{data}{the original (observed) data set.}
  \item{cont.na}{For data NOT supplied as a synthetic data object created by 
  \code{synthpop}, this gives special values for continuous variables as
   described in the documentation for the function \code{syn}.}
  \item{keys}{a vector of strings with the names of variables to be used 
    in combination to form  a quasi identifier.}   
  \item{targets}{a vector of strings with the names of variables to be used as
  	targets for the disclosure measures. Defaults to  all variables in both original
  	and synthetic data that are not in \code{keys}.}  
  \item{denom_lim}{an integer that determines the limit above which a warning 
     to check the two way relationships for potential prior disclosure information.}
  \item{exclude_ov_denom_lim}{TRUE/FALSE according to whether disclosive groups with 
    denominators > denom_lim  should be excluded from disclosure measures.} 
  \item{not.targetslev}{Vector of same length as targets giving level of each
   target to be excluded from calculating disclosure measures. Set elements for
   unaffected targets as blanks.}
  \item{print.flag}{TRUE/FALSE to print out line as disclosure for each member of
    targets is calculated.}
  \item{usetargetsNA}{A logical vector of the same length as \code{targets} that
            determines if \code{NA} values of each are to be
             considered disclosive. Defaults to \code{FAULT} for all.}
  \item{usekeysNA}{A logical vector of the same length as \code{keys} that
            determines if \code{NA} values of each key are to be
         considered disclosive. Defaults to \code{FAULT} for all keys.}
  \item{exclude.keys}{A list of same length as  \code{targets} giving the keys 
     for two way exclusions for the ith target. For details see  documentation 
     in \code{\link{disclosure}} } 
  \item{ exclude.keylevs}{A list of same length as  \code{targets} giving the levels of keys 
     for two way exclusions for the ith target. For details see  documentation 
     in \code{\link{disclosure}} } 
  \item{ exclude.targetlevs}{A list of same length as  \code{targets} giving the levels of
     target for two way exclusions for the ith target. For details see  documentation 
     in \code{\link{disclosure}} } 
  \item{ngroups_targets}{Unless set to NULL (the default) numeric target variables
     will be grouped into \code{ngroups_target} categories. If \code{ngroups_keys}
     is of length 1 all numeric  targets will be have the same number of groups. 
     Otherwise \code{ngroups_targets} needs to be a vector of the same length as 
     targets and will give the number of groups for each  target. 
     If an element of \code{ngroups_targets} is zero, no grouping will be done.}
  \item{ngroups_keys}{Unless set to NULL (the default) any numeric variable
     will be grouped into categories If \code{ngroups_keys} is of length 1 all numeric
     keys will be have the same number of groups. Otherwise \code{ngroups_keys}
     needs to be the same length as keys and will give the number of groups for each
     key. If an element of \code{ngroups_keys} is zero, no grouping will be done.}
  \item{ident.meas}{Choice of statistics to use as a measure of identity disclosure.
   Must be a selection from: \code{"repU"} or \code{"UiSiO"}. See     
        \code{\link{disclosure}} for  explanations of measures.}
  \item{attrib.meas}{Choice of statistics to use as a measure of attribute disclosure.
   Must be a selection from: \code{"DiSCO"} or \code{"DiSDiO"}. See     
        \code{\link{disclosure}} for  explanations of measures.}
  \item{thresh_1way}{A vector of two numeric values both of which meed to be exceeded
   for warnings about a level of the target that may be dominating the results. 
   The first is the count of all disclosive records, and the second is the \% of 
   all records for this level of the target. Default is c(50, 90), 
   meaning a group of 50 disclosive records for this level of the target
   where they make up over 90\% of all disclosive records.}
  \item{thresh_2way}{A vector of two numeric values both of which meed to be 
   exceeded for warnings about a level of the target that may be dominating 
   the results. The first is the count of all disclosive records for this 
   key-target combination and the second is the percantage of all disclosive 
   records for this combination. Default is c(5, 80), 
   meaning a group of more than 5 records where over 80\% of all the original 
   values with this key have this level of the target.}
  \item{digits}{number of digits to print for the disclosure measures.} 
  \item{plot}{determines if plot will be produced when the result is printed.} 
  \item{print}{logical value that determines if a summary of results is 
    to be printed.}  
  \item{compare.synorig}{a logical value  to determine if the functions
  \code{synorig.compare()} should be used to check that data sets can be 
   compared. Default set to \code{FALSE} except when the synthetic data are supplied as a    data.frame or a list when set to TRUE.} 
   \item{to.print}{Vector of items to be printed including "ident", "attrib", both or NULL}
  \item{\dots}{additional parameters}   
  \item{x}{an object of class \code{multi.disclosure}.}
}

\details{Calculates measures of identity and attribution disclosure from the keys 
  specified in \code{keys} with the function \code{\link{disclosure}}. For attribute
  disclosure a table with one line for each target  can be printed or plotted.
  Details are in help file for \code{\link{disclosure}}. 
  
}
\value{An object of class \code{multi.disclosure} which is a list with the following 
  components: 
  \item{attrib.table}{a table with the selected attribute disclosure measure
       (\code{attrib.meas}) for synthetic data and corresponding measure for the 
       original data "CAPd" if (\code{attrib.meas}) is  "DCAP", and "DiO" for others.}
  \item{attrib.plot}{plot of attrib.table with labels indicating where large
           denominators suggest checking.}
  \item{keys}{see above.}
  \item{ident.orig}{value of identity disclosure \code{UiO} from the original data, see
       help file  for \code{\link{disclosure}}.} 
  \item{ident.syn}{value of identity disclosure \code{ident.meas} from the synthetic 
  data, see help file  for \code{\link{disclosure}}.} 
  \item{Norig}{Number of records in data.}
  \item{denom_lim}{see above.}
  \item{exclude_ov_denom_lim}{see above.}
  \item{digits}{see above.} 
  \item{usetargetsNA}{see above.} 
  \item{usekeysNA}{see above.} 
  \item{ident.meas}{see above.} 
  \item{attrib.meas}{see above.} 
  \item{m}{see above.} 
  \item{plot}{see above.} 
  \item{output.list}{A named list with a component for each target 
     where each component is the output from the function 
     \code{\link{disclosure}} for that target. This
     allows check_1way and check_2way to be examined for each target.}
  \item{call}{R call used to create the object}
   }
      
\references{
to follow link to vignette
}

\seealso{
 \code{\link{disclosure}}
}

\examples{
ods <- SD2011[, c("sex", "age", "edu", "marital", "region", "income")]
s1 <- syn(ods)

### synthetic data provided as a 'data.frame' object
t1 <- multi.disclosure(s1$syn, ods,
keys = c("sex", "age", "edu"))

### synthetic data provided as a 'synds' object  
t1 <- multi.disclosure(s1, ods, 
keys = c("sex", "age", "edu"))

}
