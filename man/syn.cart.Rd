\name{syn.ctree, syn.cart}
\alias{syn.cart}
\alias{syn.ctree}
\title{Synthesis with classification and regression trees (CART)}
\description{
  Generates univariate synthetic data using classification and
  regression trees (without or with bootstrap).
}
\usage{
syn.ctree(y, x, xp, smoothing = "", proper = FALSE, 
          minbucket = 5, mincriterion = 0.9, ...)
syn.cart(y, x, xp, smoothing = "", proper = FALSE, 
         minbucket = 5, cp = 1e-08, ...)  
}
\arguments{
  \item{y}{an original data vector of length \code{n}.}
  \item{x}{a matrix (\code{n} x \code{p}) of original covariates.}
  \item{xp}{a matrix (\code{k} x \code{p}) of synthesised covariates.}
  \item{smoothing}{smoothing method for numeric variable. See 
    \code{\link{syn.smooth}}.}
  \item{proper}{for proper synthesis (\code{proper = TRUE}) a CART 
    model is fitted to a bootstrapped sample of the original data.}  
  \item{minbucket}{the minimum number of observations in
    any terminal node. See \code{\link[rpart]{rpart.control}} and
    \code{\link[partykit]{ctree_control}} for details.}
  \item{cp}{complexity parameter. Any split that does not
    decrease the overall lack of fit by a factor of cp is not
    attempted. Small values of \code{cp} will grow large trees. 
    See \code{\link[rpart]{rpart.control}} for details.}
  \item{mincriterion}{\code{1 - p-value} of the test that must be 
    exceeded for a split to be retained. Small values of 
    \code{mincriterion} will grow large trees. 
    See \code{\link[partykit]{ctree_control}} for details.}
  \item{\dots}{additional parameters passed to 
    \code{\link[partykit]{ctree_control}} for \code{syn.ctree} and 
    \code{\link[rpart]{rpart.control}} for \code{syn.cart}.}
}
\details{
  The procedure for synthesis by a CART model is 
  as follows: \enumerate{ \item Fit a classification 
  or regression tree by binary recursive partitioning. 
  \item For each \code{xp} find the terminal node. \item Randomly 
  draw a donor from the members of the node and take the observed 
  value of \code{y} from that draw as the synthetic value.}
  
  \code{syn.ctree} uses \code{\link{ctree}} function from the
  \pkg{partykit} package and \code{syn.cart} uses \code{\link{rpart}}
  function from the \pkg{rpart} package. They differ, among others, 
  in a selection of a splitting variable and a stopping rule for the 
  splitting process.
  
  A Guassian kernel smoothing can be applied to continuous variables 
  by setting smoothing parameter to \code{"density"}. It is recommended 
  as a tool to decrease the disclosure risk. Increasing \code{minbucket}
  is another means of data protection. 
  
  CART models were suggested for generation of synthetic data by 
  Reiter (2005) and then evaluated by Drechsler and Reiter (2011).
}
\value{A list with two components:
  \item{res}{a vector of length \code{k} with synthetic values of \code{y}.}
  \item{fit}{the fitted model which is an object of class \code{rpart.object} 
    or \code{ctree.object} that can be printed or plotted.}
}
\references{
  Reiter, J.P. (2005). Using CART to generate partially synthetic, 
  public use microdata. \emph{Journal of Official Statistics}, 
  \bold{21}(3), 441--462.
  
  Drechsler, J. and Reiter, J.P. (2011). An empirical evaluation of easily 
  implemented, nonparametric methods for generating synthetic datasets.
  \emph{Computational Statistics and Data Analysis}, \bold{55}(12), 3232--3243.
}
\seealso{
  \code{\link{syn}}, \code{\link{syn.survctree}},
  \code{\link[rpart]{rpart}}, \code{\link[partykit]{ctree}},
  \code{\link{syn.smooth}}
}
\keyword{datagen}
