\name{compare_canalysis}
\alias{compare_canalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compare the result of two clustering analysis }
\description{
  This function cross-compares together two cluster results from a single analysis, or two different ones. If we only consider a single analysis, then we perform cross\_comparisons of top ranking cluster results according to their BIC score. But we can also compare the result of two complete clustering analysis conducted under different settings, e.g. we might want to observe the effect of two different data preparation on the clustering. In that case, we query the function for a particular model to cross-compare.
  }
\usage{
compare_canalysis(analysis1, analysis2 = NULL, query = NULL, filePrefix = NULL, bbox_threshold = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{analysis1}{ ~~Describe \code{analysis1} here~~ }
  \item{analysis2}{ ~~Describe \code{analysis2} here~~ }
  \item{query}{ ~~Describe \code{query} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{bbox_threshold}{ ~~Describe \code{bbox_threshold} here~~ }
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

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
