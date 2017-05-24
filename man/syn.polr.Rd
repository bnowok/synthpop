\name{syn.polr}
\alias{syn.polr}
%\alias{syn.polr.proper}
\title{Synthesis by ordered polytomous regression}
\description{
  Generates a synthetic categorical variable using
  ordered polytomous regression (without or with bootstrap).
}
\usage{
syn.polr(y, x, xp, proper = FALSE, maxit = 1000, trace = FALSE, 
  MaxNWts = 10000, ...)
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{proper}{for proper synthesis (\code{proper = TRUE}) 
    a model is fitted to a bootstrapped sample of the original data.} 
  \item{maxit}{the maximum number of iterations for \code{\link[nnet]{nnet}}.}
  \item{trace}{switch for tracing optimization for \code{\link[nnet]{nnet}}.}
  \item{MaxNWts}{the maximum allowable number of weights for \code{\link[nnet]{nnet}}.}
  \item{\dots}{additional parameters passed to \code{\link[stats]{optim}} 
    or \code{\link[nnet]{nnet}}.}
}
\details{
  Generates synthetic ordered categorical variables by the
  proportional odds logistic regression (polr) model. The
  function repeatedly applies logistic regression on the
  successive splits. The model is also known as the
  cumulative link model.

  The algorithm of \code{syn.polr} uses the
  function \code{\link[MASS]{polr}} from the \pkg{MASS} package.

  In order to avoid bias due to perfect prediction, the
  data are augmented by the method of White, Daniel and Royston (2010).

  In case the call to \code{\link[MASS]{polr}} fails, 
  usually because the data are very sparse, 
  \code{\link[nnet]{multinom}} function is used instead.
}
\value{A vector of length \code{k} with synthetic values of \code{y}.}
\references{
  White, I.R., Daniel, R. and Royston, P. (2010). Avoiding bias
  due to perfect prediction in multiple imputation of
  incomplete categorical variables. \emph{Computational
  Statistics and Data Analysis}, \bold{54}, 2267--2275.
}
\seealso{
  \code{\link{syn}},\code{\link{syn.polyreg}}
  \code{\link[nnet]{multinom}},
  \code{\link[MASS]{polr}}
}
\keyword{datagen}

