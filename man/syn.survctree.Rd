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
    any terminal node. See \code{\link{ctree_control}}
    for details.}
  \item{\dots}{additional parameters passed to \code{\link{ctree}}.}
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
  
  NOTE that when the function is called by setting elements of method in 
  \code{syn()} to \code{"survctree"}, the parameter \code{minbucket}
  can be supplied to \code{syn()} as e.g. \code{survctree.minbucket}.
}
\value{
  A list with the following components:
  \item{syn.time}{a vector of length \code{k} with synthetic time values.}
  \item{syn.event}{a vector of length \code{k} with synthetic event indicator values.}
}
%references{...}
\seealso{
  \code{\link{syn}}, \code{\link{syn.ctree}}
}
\keyword{datagen}
