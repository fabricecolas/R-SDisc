\name{stability_assessment}
\alias{stability_assessment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stability_assessment(data = na.omit(as.data.frame(df)), analysis_outcomes = NULL, sumscore_groups = NULL, stats_fun = NULL, fun_transform = list(center = transform_AVG, scale = transform_SIGMA), G_set = 1:9, query_model = "VVI,5", modelNames = c("EII", "VII", "EEI", "EVI", "VEI", "VVI", "EEE", "EEV", "VEV", "VVV"), noise_vector = (1/2^(1:10)), K = 10, filePrefix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{analysis_outcomes}{ ~~Describe \code{analysis_outcomes} here~~ }
  \item{sumscore_groups}{ ~~Describe \code{sumscore_groups} here~~ }
  \item{stats_fun}{ ~~Describe \code{stats_fun} here~~ }
  \item{fun_transform}{ ~~Describe \code{fun_transform} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{query_model}{ ~~Describe \code{query_model} here~~ }
  \item{modelNames}{ ~~Describe \code{modelNames} here~~ }
  \item{noise_vector}{ ~~Describe \code{noise_vector} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
}
\details{
  
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
