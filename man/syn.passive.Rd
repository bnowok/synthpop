\name{syn.passive}
\alias{syn.passive}
\title{Passive synthesis}
\description{Derives a new variable according to a specified 
  function of synthesised data.}
\usage{syn.passive(data, func)}
\arguments{
  \item{data}{a data frame with synthesised data.} 
  \item{func}{a \code{formula} specifying transformations on data. 
    It is specified as a string starting with \code{~}.}
}
\details{
  Any function of the synthesised data can be specified. Note that several operators such as 
  \code{+}, \code{-}, \code{*} and \code{^} have different meanings in \code{formula} syntax. Use 
  the identity function \code{I()} if they should be interpreted as arithmetic operators, 
  e.g. \code{"~I(age^2)"}.
  }
\value{A vector including the result of applying the \code{formula}.}
\references{
  Van Buuren, S. and Groothuis-Oudshoorn, K. (2011).
  \code{mice}: Multivariate Imputation by Chained Equations
  in \code{R}. \emph{Journal of Statistical Software},
  \bold{45}(3), 1-67.
  \url{http://www.jstatsoft.org/v45/i03/}
}
\author{Stef van Buuren, Karin Groothuis-Oudshoorn, 2000}
\seealso{
  \code{\link{syn}}
}
\keyword{datagen}
