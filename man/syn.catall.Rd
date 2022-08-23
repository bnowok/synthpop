\name{syn.catall}
\alias{syn.catall}
\title{Synthesis of a group of categorical variables from a saturated model}
\description{
  A saturated model is fitted to a table produced by cross-tabulating
  all the variables.}
\usage{
syn.catall(x, k, proper = FALSE, priorn = 1, structzero = NULL, 
           maxtable = 1e8, epsilon = 0, rand = TRUE,  ...)}
\arguments{
  \item{x}{a data frame (\code{n} x \code{p}) of the set of original variables.}
  \item{k}{a number of rows in each synthetic data set - defaults to \code{n}.}
  \item{proper}{if \code{proper = TRUE} \code{x} is replaced with a bootstrap 
    sample before synthesis, thus effectively sampling from the posterior 
    distribution of the model, given the data.}
  \item{priorn}{the sum of the parameters of the Dirichelet prior which can be 
    thought of as a pseudo-count giving the number of observations that inform 
    prior knowledge about the parameters.}
  \item{structzero}{a named list of lists that defines which cells in the table 
    are structural zeros and will remain as zeros in the synthetic data, by 
    leaving their prior as zeros. Each element of the \code{structzero} list is 
    a list that describes a set of cells in the table defined by a combination 
    of two or more variables and a name of each such element must consist of 
    those variable names seperated by an underscore, e.g. \code{sex_edu}. 
    The length of each such element is determined by the number of variables and 
    each component gives the variable levels (numeric or labels) that define the 
    structural zero cells (see an example below).}
  \item{maxtable}{a number of cells in the cross-tabulation of all the variables
    that will trigger a severe warning.}
  \item{epsilon}{measures scale of laplace noise to be added under differential 
    privacy (DP)}
  \item{rand}{for DP versions determines if multinomial noise is to be added to 
    DP counts. If it is set to false the DP adjusted counts are simply rounded 
    to a whole number in a manner that preserves the desired sample size (k).}
  \item{\dots}{additional parameters.}
}
\details{
When used in \code{syn} function the group of categorical variables 
with \code{method = "catall"} must all be together at the start of the 
\code{visit.sequence}. Subsequent variables in \code{visit.sequence} are then 
synthesised conditional on the synthesised values of the grouped variables.
A saturated model is fitted to a table produced by cross-tabulating all the 
variables. Prior probabilities for the proportions in each cell of the table
are specified from the parameters of a Dirichlet distribution with the same 
parameter for every cell in the table that is not a structural zero (see above).
The sum of these parameters is \code{priorn} so that each one is \eqn{priorn/N}
where \eqn{N} is the number of cells in the table that are not structural zeros.
The default \code{priorn = 1} can be thought of as equivalent to the knowledge 
that \code{1} observation would be equally likely to be in any cell that is not 
a structural zero. The posterior expectation, given the observed counts,
for the probability of being in a cell with observed count \eqn{n_i}
is thus \eqn{(n_i + priorn/N) / (N + priorn)}. The synthetic data are generated 
from a multinomial distribution with parameters given by these probabilities.

Unlike \code{syn.satcat}, which fits saturated conditional models,
the synthesised data can include any combination of variables, except 
those defined by the combinations of variables in \code{structzero}.

NOTE that when the function is called by setting elements of method in 
\code{syn()} to \code{"catall"}, the parameters \code{priorn}, \code{structzero}, 
\code{maxtable}, \code{epsilon}, and \code{rand} must be supplied to \code{syn} 
as e.g. \code{catall.priorn}.
}
\value{
A list with two components:
  \item{res}{a data frame of dimension \code{k x p} containing the synthesised data.}
  \item{fit}{the cross-tabulation of all the original variables used.}
}
\examples{
ods <- SD2011[, c(1, 4, 5, 6, 2, 10, 11)]
table(ods[, c("placesize", "region")])

# Each `placesize_region` sublist: 
# for each relevant level of `placesize` defined in the first element, 
# the second element defines regions (variable `region`) that do not 
# have places of that size.

struct.zero <- list(
  placesize_region = list(placesize = "URBAN 500,000 AND OVER", 
                          region = c(2, 4, 5, 8:13, 16)),
  placesize_region = list(placesize = "URBAN 200,000-500,000", 
                          region = c(3, 4, 10:11, 13)),
  placesize_region = list(placesize = "URBAN 20,000-100,000", 
                          region = c(1, 3, 5, 6, 8, 9, 14:15)))

syncatall <- syn(ods, method = c(rep("catall", 4), "ctree", "normrank", "ctree"),
                 catall.priorn = 2, catall.structzero = struct.zero)
}
\keyword{datagen}
