\name{retrieve_model}
\alias{retrieve_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
retrieve_model(canalysis, query)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{canalysis}{ ~~Describe \code{canalysis} here~~ }
  \item{query}{ ~~Describe \code{query} here~~ }
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
function(canalysis,query){
   # init
   params       <- list()
   # retrieve the model corresponding to the query
   mbc_params   <- strsplit(query,",")[[1]]
   params       <- list(model=mbc_params[1],G=as.numeric(mbc_params[2]))
   i            <- 1 
   while(!(canalysis[["out"]][[i]]$G == params$G & canalysis[["out"]][[i]]$modelName == params$model & i <= length(canalysis[["out"]])))
      i         <- i+1
   params[["model"]]     <- canalysis[["out"]][[i]]
   params[["labelling"]] <- map(params[["model"]]$z,warn=FALSE)
   return(params)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
