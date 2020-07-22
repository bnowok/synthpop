\name{syn.polyreg}
\alias{syn.polyreg}
%\alias{syn.polyreg.proper}
\title{Synthesis by unordered polytomous regression}
\description{
  Generates a synthetic categorical variable using
  unordered polytomous regression (without or with bootstrap).
}
\usage{
syn.polyreg(y, x, xp, proper = FALSE, maxit = 1000, trace = FALSE, 
            MaxNWts = 10000, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{proper}{for proper synthesis (\code{proper = TRUE}) 
    a multinomial model is fitted to a bootstrapped sample of the original data.} 
  \item{maxit}{the maximum number of iterations for \code{\link[nnet]{nnet}}.}
  \item{trace}{switch for tracing optimization for \code{\link[nnet]{nnet}}.}
  \item{MaxNWts}{the maximum allowable number of weights for \code{\link[nnet]{nnet}}.}
  \item{\dots}{additional parameters passed to \code{\link[nnet]{nnet}}.}
}
\details{
  Generates synthetic categorical variables by the polytomous 
  regression model. The method consists of the following steps: \enumerate{
  \item Fit categorical response as a multinomial model.
  \item Compute predicted categories. \item Add appropriate
  noise to predictions.}
  
  The algorithm of \code{syn.polyreg} uses the function 
  \code{\link[nnet]{multinom}} from the \pkg{nnet} package. Any numerical 
  variables are scaled to cover the range (0,1) before fitting. Warnings 
  are printed if the algorithm fails to converge in \code{maxit} iterations
  and also if the synthesised data has only one category. The latter may occur 
  if the variable being synthesised is sparse so that the algorithm fails to 
  iterate.

  In order to avoid bias due to perfect prediction, the
  data are augmented by the method of White, Daniel and Royston (2010).   
  
  NOTE that when the function is called by setting elements of method in \code{syn()}
  to \code{"polyreg"}, the parameters \code{maxit}, \code{trace} and \code{MaxNWts}
  can be supplied to \code{syn()} as e.g. \code{polyreg.maxit}.
}
\value{
A list with two components:
  \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
  \item{fit}{a summary of the model fitted to the observed data and used
    to produce synthetic values.}
}
\references{
  White, I.R., Daniel, R. and Royston, P. (2010). Avoiding bias
  due to perfect prediction in multiple imputation of
  incomplete categorical variables. \emph{Computational
  Statistics and Data Analysis}, \bold{54}, 2267--2275.}
\seealso{
  \code{\link{syn}}, \code{\link{syn.polr}},
  \code{\link[nnet]{multinom}}, \code{\link[MASS]{polr}}
}
\keyword{datagen}

