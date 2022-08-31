\name{syn.survctree}
\alias{syn.survctree}
%\alias{syn.survctree.proper}
\title{Synthesis of survival time by classification and regression trees (CART)}
\description{
  Generates synthetic event indicator and time to event data using 
  classification and regression trees (without or with bootstrap).
}
\usage{
syn.survctree(y, yevent, x, xp, proper = FALSE, minbucket = 5, ...)
}
\arguments{
  \item{y}{a vector of length \code{n} with original time data.}
  \item{yevent}{a vector of length \code{n} with original event indicator data.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{proper}{for proper synthesis (\code{proper = TRUE}) a CART 
    model is fitted to a bootstrapped sample of the original data.}  
  \item{minbucket}{the minimum number of observations in
    any terminal node. See \code{\link[party]{ctree_control}}
    for details.}
  \item{\dots}{additional parameters passed to \code{\link[party]{ctree}}.}
}
\details{
  The procedure for synthesis by a CART model is as follows: 
  \enumerate{\item Fit a tree-structured survival model by 
  binary recursive partitioning (the terminal nodes include 
  Kaplan-Meier estimates of the survival time). 
  \item For each \code{xp} find the terminal node. \item Randomly 
  draw a donor from the members of the node and take the observed 
  value of \code{yevent} and \code{y} from that draw as the 
  synthetic values.}    
  
  The function is used in \code{syn()} to generate survival times 
  by setting elements of method in \code{syn()} to \code{"survctree"}.
  Additional parameters related to \code{\link[party]{ctree}} function, 
  e.g. \code{minbucket} can be supplied to \code{syn()} as 
  \code{survctree.minbucket}.
  
  Where the survival variable is censored this information must be supplied 
  to \code{syn()} as a named list (event) that gives the name of the variable 
  for each event indicator. Event variables can be a numeric variable with 
  values 1/0 (1 = event), TRUE/FALSE (TRUE = event) or a factor with 2 levels 
  (level 2 = event). The event variable(s) will be synthesised along with the 
  survival time(s).
}

\value{
  A list with the following components:
    \item{syn.time}{a vector of length \code{k} with synthetic time values.}
    \item{syn.event}{a vector of length \code{k} with synthetic event indicator
      values.}
    \item{fit}{the fitted model which is an item of class \code{ctree.object}.}
}

\seealso{
  \code{\link{syn}}, \code{\link{syn.ctree}}
}

\examples{
### This example uses the data set 'mgus2' from the survival package.
### It has a follow-up time variable 'futime' and an event indicator 'death'.
library(survival)

### first exclude the 'id' variable and run a dummy synthesis to get 
### a method vector
ods <- mgus2[-1]
s0 <- syn(ods)

### create new method vector including 'survctree' for 'futime' and create 
### an event list for it; the names of the list element must correspond to 
### the name of the follow-up variable for which the event indicator
### need to be specified.
meth <- s0$method
meth[names(meth) == "futime"] <- "survctree"
evlist <- list(futime = "death")

s1 <- syn(ods, method = meth, event = evlist)

### evaluate outputs
## compare selected variables
compare(s1, ods, vars = c("futime", "death", "sex", "creat"))

## compare original and synthetic follow up time by an event indicator
multi.compare(s1, ods, var = "futime", by = "death")

## compare survival curves for original and synthetic data
par(mfrow = c(2,1))
plot(survfit(Surv(futime, death) ~ sex, data = ods), 
     col = 1:2, xlim = c(0,450), main = "Original data")
legend("topright", levels(ods$sex), col = 1:2, lwd = 1, bty = "n")
plot(survfit(Surv(futime, death) ~ sex, data = s1$syn), 
     col = 1:2, xlim = c(0,450), main = "Synthetic data")
}

\keyword{datagen}
