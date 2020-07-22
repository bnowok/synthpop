\name{summary.fit.synds}
\alias{summary.fit.synds}
\alias{print.summary.fit.synds}
\title{Inference from synthetic data}
\description{Combines the results of models fitted to each of the \code{m} 
  synthetic data sets.}
\usage{
\method{summary}{fit.synds}(object, population.inference = FALSE, msel = NULL,
  real.varcov = NULL, \dots)

\method{print}{summary.fit.synds}(x, \dots)
}

\arguments{
 \item{object}{an object of class \code{fit.synds} created by fitting a model to 
   synthesised data set using function \code{\link{glm.synds}}, 
   \code{\link{lm.synds}},\code{\link{multinom.synds}} or \code{\link{polr.synds}}.}
 \item{population.inference}{a logical value indicating whether inference 
   should be made to population quantities. If \code{FALSE} inference is made
   to the results that would be expected from an analysis of the original data. 
   This option should be selected if the synthetic data are being used for 
   exploratory analysis, but the final published results will be obtained by 
   running code on the original confidential data. If \code{population.inference = TRUE} 
   results would allow population inference to be made from the synthetic data. 
   In both cases the inference will depend on the synthesising model being correct,
   but this can be checked by running the same analysis on the real data,
   see \code{\link{compare.fit.synds}}.}
 \item{msel}{index or indices of the synthetic datasets (\code{1}, \code{...}, \code{m}), 
   for which summaries of fitted models are to be produced. If \code{NULL} (default) 
   only the summary of combined estimates is produced.} 
 \item{real.varcov}{the estimated variance-covariance matrix of the fit of the model 
   to the original data. This parameter is used in the function \code{compare.fit.synds} 
   which has the original data as one of its parameters.}
 \item{\dots}{additional parameters.}
 \item{x}{an object of class \code{summary.fit.synds}.}
}
\details{The mean of the estimates from each of the m synthetic data sets 
yields asymptotically unbiased estimates of the coefficients if the observed 
data conform to the distribution used for synthesis. The standard errors are 
estimated differently depending whether inference is made for the results that 
we would expect to obtain from the observed data or for the parameters of the 
population that we assume the observed data are sampled from. The standard 
errors also differ according to whether synthetic data were produced using 
simple or proper synthesis (for details see Raab et al. (2017)).}

\value{
  An object of class \code{summary.fit.synds} which is a list with the following 
  components:
  \item{call}{the original call to \code{glm.synds} or \code{lm.synds}.}
  \item{proper}{a logical value indicating whether synthetic data were 
    generated using proper synthesis.}
  \item{population.inference}{a logical value indicating whether inference 
    is made to population coefficients or to the results that would be expected 
    from an analysis of the original data (see above).} 
  \item{incomplete}{a logical value indicating whether any of the variables in 
    the model were not synthesised. It is derived in the synthpop implementation 
    of the fitting functions (\code{\link{lm.synds}}, \code{\link{glm.synds}},
    \code{\link{multinom.synds}} and \code{\link{polr.synds}}) and saved with the 
    fitted object. When \code{TRUE} inference with \code{population.inference = TRUE}
    uses the method proposed by Reiter (2003) for what he terms partially synthetic 
    data. This method requires multiple syntheses (\code{m > 1}). If \code{m = 1},
    \code{incomplete = TRUE} and \code{population.inference = TRUE} the results will
    be still calculated and returned with warning. This will usually give standard 
    errors that are larger than they should be.} 
  \item{fitting.function}{function used to fit the model.} 
  \item{m}{the number of synthetic versions of the original (observed) data.} 
  \item{coefficients}{a matrix with combined estimates. If inference is required
    to the results that would be obtained from an analysis of the original data, 
    (\code{population.inference = FALSE}) the coefficients are given by 
    \code{xpct(Beta)}, the standard errors by \code{xpct(se.Beta)} and the 
    corresponding Z-statistic by \code{xpct(Z)}. If the synthetic data are to 
    be used to make inferences to population quantities 
    (\code{population.inference = TRUE}), the coefficients are given by 
    \code{Beta.syn}, their standard errors by \code{se.Beta.syn} and the 
    Z-statistic by  \code{Z.syn} (see vignette on inference for more details).}
  \item{n}{a number of cases in the original data.} 
  \item{k}{the number of cases in the synthesised data. Note that if \code{k} 
    and \code{n} are not equal and \code{population.inference = FALSE} (the default), 
    then the standard errors produced will estimate what would be expected by 
    an analysis of the original data set of size \code{n}.}
  \item{analyses}{\code{summary.glm} or \code{summary.lm} object respectively
    or a list of \code{m} such objects.} 
  \item{msel}{index or indices of synthetic data copies for which summaries 
    of fitted models are produced. If \code{NULL} only a summary of combined 
    estimates is produced.}
}

\references{
  Nowok, B., Raab, G.M and Dibben, C. (2016). synthpop: Bespoke
  creation of synthetic data in R. \emph{Journal of Statistical Software},
  \bold{74}(11), 1-26. \doi{10.18637/jss.v074.i11}.

  Raab, G.M., Nowok, B. and Dibben, C. (2017). Practical data synthesis for 
  large samples. \emph{Journal of Privacy and Confidentiality}, \bold{7}(3), 
  67-97. Available at: 
  \url{https://journalprivacyconfidentiality.org/index.php/jpc/article/view/407}

  Reiter, J.P. (2003) Inference for partially synthetic, public use microdata 
  sets. \emph{Survey Methodology}, \bold{29}, 181-188.
}
\seealso{
\code{\link{compare.fit.synds}}, \code{\link{summary}}, \code{\link{print}} 
}
\examples{
ods <- SD2011[1:1000,c("sex","age","edu","ls","smoke")]
  
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
