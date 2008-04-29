\name{get_coloring_scheme}
\alias{get_coloring_scheme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_coloring_scheme(data, canalysis_variables = parkinson_canalysis_variables)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{canalysis_variables}{ ~~Describe \code{canalysis_variables} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
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
function(data, canalysis_variables = parkinson_canalysis_variables){
   data <- data.matrix(data)
   # MAKE A DENDROGRAM FROM THESE AVG PATTERN FOR THE DENDROGRAM ORDERING
   dendro_patterns      <- hclust(dist(data))
   # SELECT DIVERGING COLORS FROM RColorBrewer
   cluster_colors       <- brewer.pal(nrow(data),"Set1")
   # RETURN A DATA STRUCTURE WITH CLUSTERS RIGHTLY ORDERED AND THEIR ASSOCIATED COLOR
   return(list(cluster_order=dendro_patterns$order,cluster_color=cluster_colors))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
