\name{summary.fit.synds}
\alias{summary.fit.synds}
\alias{print.summary.fit.synds}
\title{Inference from synthetic data}
\description{Combines the results of models fitted to each of the \code{m} 
  synthetic data sets.}
\usage{
\method{summary}{fit.synds}(object, population.inference = FALSE, msel = NULL,
partly = FALSE, \dots)

\method{print}{summary.fit.synds}(x, \dots)
}

\arguments{
  \item{object}{an object of class \code{fit.synds} created by fitting a model to 
   synthesised data set using function \code{\link{glm.synds}} or
   \code{\link{lm.synds}}.}
  \item{population.inference}{a logical value indicating whether inference 
   should be made to population quantities. If \code{FALSE} inference is made
   to original data quantities.}
  \item{msel}{index or indices of synthetic data copies for which summaries 
    of fitted models are to be produced. If \code{NULL} (default) a summary 
    of combined estimates is produced.}
  \item{partly}{a logical value indicating whether data are partly synthesised.}
  \item{\dots}{additional parameters.}
  \item{x}{an object of class \code{summary.fit.synds}.}
}
\details{The mean of the estimates from each of the m synthetic data sets 
yields asymptotically unbiased estimates of the coefficients if the observed 
data conform to the distribution used for synthesis. The standard errors are 
estimated differently depending whether inference is made for the results that 
would be obtained from the observed data or for the parameters of the population 
that we assume the observed data are sampled from. The standard errors also differ 
according to whether synthetic data were produced using simple or proper 
synthesis (for details see Raab et al. (submitted 2014)).}

\value{
  An object of class \code{summary.fit.synds} which is a list with the following 
  components:
  \item{call}{the original call to \code{glm.synds} or \code{lm.synds}.}
  \item{proper}{a logical value indicating whether synthetic data were 
    generated using proper synthesis.}
  \item{population.inference}{a logical value indicating whether inference to 
    population coefficients or to coefficients of the actual (observed) data is made.} 
  \item{fitting.function}{function used to fit the model.} 
  \item{m}{the number of synthetic versions of the original (observed) data.} 
  \item{coefficients}{a matrix with combined estimates. It includes point 
   estimates of coefficients (\code{B.syn}), their standard errors 
   (\code{se(B.syn)}) and Z scores (\code{Z.syn}) for population 
   and observed data quantities respectively. For inference to original data 
   quantities it contains in addition estimates of the actual standard errors 
   based on synthetic data (\code{se(Beta).syn}) and standard errors of 
   Z scores (\code{se(Z.syn)}).} 
  \item{n}{a number of cases in the original data.} 
  \item{k}{a number of cases in the synthesised data.}
  \item{analyses}{\code{summary.glm} or \code{summary.lm} object respectively
    or a list of \code{m} such objects.} 
  \item{msel}{index or indices of synthetic data copies for which summaries 
    of fitted models are produced. If \code{NULL} a summary of combined 
    estimates is produced.}
}

\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.

  Raab, G.M., Nowok, B. and Dibben, C. (resubmitted 2016). 
  Practical data synthesis for large samples. \url{https://arxiv.org/pdf/1409.0217v7}
}
\seealso{
\code{\link{summary}},\code{\link{print}}
}
\examples{
ods <- SD2011[1:2000,c("sex","age","edu","ls","smoke")]
  
### simple synthesis
s1 <- syn(ods, m = 5)
f1 <- glm.synds(smoke ~ sex + age + edu + ls, data = s1, family = "binomial")
summary(f1)
summary(f1, population.inference = TRUE)
  
### proper synthesis
s2 <- syn(ods, m = 5, method = "parametric", proper = TRUE)
f2 <- glm.synds(smoke ~ sex + age + edu + ls, data = s2, family = "binomial")
summary(f2)
summary(f2, population.inference = TRUE)
}
