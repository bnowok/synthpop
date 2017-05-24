\name{SD2011}
\alias{SD2011}
\docType{data}
\title{
  Social Diagnosis 2011 - Objective and Subjective Quality of Life in Poland
}
\description{
  Sample of 5,000 individuals from the Social Diagnosis 2011 survey;
  selected variables only.
}
\usage{SD2011}
\format{A data frame with 5,000 observations on the following 35 variables:
\describe{
\item{sex}{Sex}
\item{age}{Age of person, 2011}
\item{agegr}{Age group, 2011}
\item{placesize}{Category of the place of residence}
\item{region}{Region (voivodeship)}
\item{edu}{Highest educational qualification, 2011}
\item{eduspec}{Discipline of completed qualification}
\item{socprof}{Socio-economic status, 2011}
\item{unempdur}{Total duration of unemployment in the last 2 years (in months)}
\item{income}{Personal monthly net income}
\item{marital}{Marital status}
\item{mmarr}{Month of marriage}
\item{ymarr}{Year of marriage}
\item{msepdiv}{Month of separation/divorce}
\item{ysepdiv}{Year of separation/divorce}
\item{ls}{Perception of life as a whole}
\item{depress}{Depression symptoms indicator}
\item{trust}{View on interpersonal trust}
\item{trustfam}{Trust in own family members}
\item{trustneigh}{Trust in neighbours}
\item{sport}{Active engagement in some form of sport or exercise}
\item{nofriend}{Number of friends}
\item{smoke}{Smoking cigarettes}
\item{nociga}{Number of cigarettes smoked per day}
\item{alcabuse}{Drinking too much alcohol}
\item{alcsol}{Starting to use alcohol to cope with troubles}
\item{workab}{Working abroad in 2007-2011}
\item{wkabdur}{Total time spent on working abroad}
\item{wkabint}{Plans to go abroad to work in the next two years}
\item{wkabintdur}{Intended duration of working abroad}
\item{emcc}{Intended destination country}
\item{englang}{Knowledge of English language}
\item{height}{Height of person}
\item{weight}{Weight of person}
\item{bmi}{Body mass index}
}}
\note{Please note that the original variable names have been changed to make them 
  more self-explanatory. Some variable labels have been adjusted as well.}
\source{
  Council for Social Monitoring. Social Diagnosis 2000-2011: integrated database.
  \url{http://www.diagnoza.com/index-en.html} [downloaded on 13/12/2013]
}
\references{
  Czapinski J. and Panek T. (Eds.) (2011). Social Diagnosis 2011. Objective and 
  Subjective Quality of Life in Poland - full report. Contemporary Economics, 
  Volume 5, Issue 3 (special issue) \url{http://ce.vizja.pl/en/issues/volume/5/issue/3#art254} 
}
\examples{
  spineplot(englang ~ agegr, data = SD2011, xlab = "Age group", ylab = "Knowledge of English")
  boxplot(income ~ sex, data = SD2011[SD2011$income != -8,])
}

\keyword{datasets}
