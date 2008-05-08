\name{analysis}
\alias{analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
analysis(data, config, calgo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{config}{ ~~Describe \code{config} here~~ }
  \item{calgo}{ ~~Describe \code{calgo} here~~ }
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
function(data, config, calgo){
   # REDIRECT GRAPHICS -> PS
   prefix       <- attr(config,"prefix")
   # ------------
   postscript(file=paste2(prefix,"_Cluster_Analysis.ps"))
      # | VISUALIZE DISTRIB AND CORRELATION ON ORIGINAL AND TRANSFORMED DATA
      variable_graphic_report(data, which_data = "original", prefix = prefix)
      variable_graphic_report(data, prefix = prefix)
      # | CLUSTERING ANALYSIS 
      cresult           <- calgo(data, config)
      # | VISUAL CHARACTERISTICS OF THE CLUSTERS
      cresult           <- cresult_graphic_characteristics(cresult,canalysis_data = data,canalysis_config = config)
   graphics.off()
   # ------------
   # CLUSTER RESULT CHARACTERISTIZATION AND EVALUATION
   cresult              <- cresult_numeric_characteristics(cresult, canalysis_data, canalysis_config)
   # COMPARE BEST CLUSTERS  2-BY-2 AND SAVE OUTPUT INTO CSV FILES
   compare_canalysis(cresult, config=config)
   return(cresult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
