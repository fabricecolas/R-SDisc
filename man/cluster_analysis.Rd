\name{cluster_analysis}
\alias{cluster_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generic cluster analysis function}
\description{
We define that function to carry a generic cluster analysis by means of some clustering algorithm. In our case, we only proceed to a cluster analysis by mixture modeling, see model\_based\_clustering for further details. In order to summary and characterize a cluster result, we may provide a list of functions as parameters in "stats\_fun", which will derive additional statistics on the clustered data.
}
\usage{
cluster_analysis(data, filePrefix = filePrefix, G_set = 1:9, stats_fun = list(oddratios = function(data, class) {
    return(statistics_logodds(data, class))
}))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{stats_fun}{ ~~Describe \code{stats_fun} here~~ }
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
