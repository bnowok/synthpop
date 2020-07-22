\name{syn.nested}
\alias{syn.nested}
\title{Synthesis for a variable nested within another variable.}
\description{Synthesizes one variable (\code{y}) from another one (\code{x}) 
  when \code{y} is nested in the categories of \code{x}. A bootstrap 
  sample is created from the original values of \code{y} within each category 
  of \code{xp} (the synthesised values of the grouping variable).}
\usage{
syn.nested(y, x, xp, smoothing = "", cont.na = NA,  ...) 
}
\arguments{
  \item{y}{an original data vector of length \code{n} for the nested variable.}
  \item{x}{an original data vector of length \code{n} for the variable within 
    which \code{y} is nested.}
  \item{xp}{a vector of length \code{k} with synthetic values of \code{x}.}
  \item{smoothing}{if this has the value "density"  then  values of a numeric 
    \code{y} are smoothed after synthesis with a kernel density smoother.}  
  \item{cont.na}{when y is numeric this can be a list or a vector giving values 
    of \code{y} that indicate missing values.}
  \item{\dots}{additional parameters.}
}
\details{An example would be when \code{x} is a classification 
  of occupations and \code{y} is a more detailed sub-classification. It is 
  intended that \code{x} is a categorical (factor) variable.
  A warning will be issued if the original \code{y} is not nested within \code{x}.
  A variable synthesised by \code{syn.nested()} is automatically excluded from
  predicting later variables because it will provide no extra information,
  given its grouping variable.
  \code{syn.nested()} is also used for the final synthesis of variables in 
  \code{syn()} when the option \code{numtocat} is used to synthesise numerical 
  variables as groups.}
\value{
A list with two components:
  \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
  \item{fit}{a name of the method used for synthesis (\code{"nested"}).}
}
\keyword{datagen}
