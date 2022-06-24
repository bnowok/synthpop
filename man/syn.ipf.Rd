\name{syn.ipf}
\alias{syn.ipf}
\title{Synthesis of a group of categorical variables by iterative 
  proportional fitting}
\description{
  A fit to the table is obtained from the log-linear fit that matches the 
  numbers in the margins specified by the margin parameters.}
\usage{
syn.ipf(x, k, proper = FALSE, priorn = 1, structzero = NULL, 
        gmargins = "twoway", othmargins = NULL, tol = 1e-4,
        max.its = 5000, maxtable = 1e8, print.its = FALSE,
        epsilon = 0, rand = TRUE, ...)
}
\arguments{
  \item{x}{a data frame of the set of original data to be synthesised.}
  \item{k}{a number of rows in each synthetic data set - defaults to \code{n}.}
  \item{proper}{if \code{proper = TRUE} \code{x} is replaced with a bootstrap 
    sample before synthesis, thus effectively sampling from the posterior 
    distribution of the model, given the data.}
  \item{priorn}{the sum of the parameters of the Dirichlet prior which can be 
    thought of as a pseudo-count giving the number of observations that inform 
    prior knowledge about the parameters.}
  \item{structzero}{NOT IMPLEMENTED YET
    a named list of lists that defines which cells in the table 
    are structural zeros and will remain as zeros in the synthetic data, by 
    leaving their prior as zeros. Each element of the \code{structzero} list is 
    a list that describes a set of cells in the table defined by a combination 
    of two or more variables and a name of each such element must consist of 
    those variable names seperated by an underscore, e.g. \code{sex_edu}. 
    The length of each such element is determined by the number of variables and 
    each component gives the variable levels (numeric or labels) that define the 
    structural zero cells (see an example below).}
  \item{gmargins}{a single character to define a group of margins. At present 
    there is "oneway" and "twoway" option that creates, respectively, all 1-way 
    and 2-way margins from the table.}
  \item{othmargins}{a list of margins that will be fitted. If \code{gmargins} is
    not \code{NULL} \code{othmargins} will be added to them.}
  \item{tol}{stopping criterion for \code{\link{Ipfp}}.}
  \item{max.its}{maximum umber of iterations allowed for \code{\link{Ipfp}}.}
  \item{maxtable}{the number of cells in the cross-tabulation of all the 
    variables that will trigger a severe warning.}
  \item{print.its}{if true the iterations from \code{\link{Ipfp}} will be 
    printed on the console. Otherwise only a message as to whether the 
    iterations have converged will be given at the end of the fitting.}
  \item{epsilon}{epsilon value for overall differential privacy (DP) parameter.
    This is implemented by dividing the privacy budget equally over all the 
    margins used to fit the data.}  
  \item{rand}{when epsilon is > 0 and DP synthetic data are created this 
    determines whether the data are created by Poisson counts from the expected 
    fitted counts in the cells of the DP adjusted data.}
  \item{\dots}{additional parameters.}
}
\details{
When used in \code{syn} function the group of variables with 
\code{method = "ipf"} must all be together at the start of the visit sequence. 
This function is designed for categorical variables, but it can also be used for 
numerical variables if they are categorised by specifying them in the 
\code{numtocat} parameter of the main function \code{syn}. Subsequent variables 
in \code{visit.sequence} are then synthesised conditional on the synthesised 
values of the grouped variables. A fit to the table is obtained from the 
log-linear fit that matches the numbers in the margins specified by the margin 
parameters. Prior probabilities for the proportions in each cell of the table 
are given by a Dirichlet distribution with the same parameter for every cell 
in the table that is not a structural zero. The sum of these parameters is 
\code{priorn}. The default \code{priorn = 1} can be thought of as equivalent 
to the knowledge that \code{1} observation would be equally likely to
fall in any cell of the table. The synthetic data are generated from a multinomial
distribution with parameters given by the expected posterior probabilities for 
each cell of the table. If the maximum likelihood estimate from the log-linear 
fit to cell \eqn{c_i} is \eqn{p_i} and the table has \eqn{N} cells that are not
structural zeros then the expectation of the posterior probability
for this cell is \eqn{(p_i + priorn/N^2) / (1 + priorn / N^2)} or
equivalently \eqn{(N * p_i + priorn/N) / (N + priorn / N)}.

Unlike \code{syn.satcat}, which fits saturated models from their conditional 
distrinutions, \code{x} can include any combination of variables, including 
those not present in the original data, except those defined by \code{structzero}.

NOTE that when the function is called by setting elements of
method in \code{syn} to \code{"ipf"}, the parameters \code{priorn},
\code{structzero}, \code{gmargins}, \code{othmargins}, \code{tol}, 
\code{max.its}, \code{maxtable}, \code{print.its}, \code{epsilon}, 
and \code{rand} must be supplied to \code{syn} as e.g. \code{ipf.priorn}.
}
\value{
A list with two components:
  \item{res}{a data frame with \code{k} rows containing the synthesised data.}
  \item{fit}{a list made up of two lists: the margins fitted and the original 
    data for each margin.}
}
\examples{
ods <- SD2011[, c(1, 4, 5, 6, 2, 10, 11)]
table(ods[, c("placesize", "region")])

synipf <- syn(ods, method = c(rep("ipf", 4), "ctree", "normrank", "ctree"), 
              ipf.gmargins = "twoway", ipf.othmargins = list(c(1, 2, 3)),
              ipf.priorn = 2)
}
\keyword{datagen}
