\name{analysis}
\alias{analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate visual-aids and table-charts in a cluster analysis with different settings}
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analysis(data, analysis_outcomes = list(), sumscore_groups = list(), stats_fun = list(), G_set = 1:9, formula_matrix = NULL, fun_transform = list(identity = function(x) {
    return(data = I(x), model = NULL)
}), filePrefix = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{analysis_outcomes}{ ~~Describe \code{analysis_outcomes} here~~ }
  \item{sumscore_groups}{ ~~Describe \code{sumscore_groups} here~~ }
  \item{stats_fun}{ ~~Describe \code{stats_fun} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{formula_matrix}{ ~~Describe \code{formula_matrix} here~~ }
  \item{fun_transform}{ ~~Describe \code{fun_transform} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
