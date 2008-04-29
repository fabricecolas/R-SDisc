\name{em_mbc}
\alias{em_mbc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em_mbc(data, modelName = "VVI", G = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{modelName}{ ~~Describe \code{modelName} here~~ }
  \item{G}{ ~~Describe \code{G} here~~ }
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
function(data, modelName="VVI", G=5) {
   if(G >= 2){
      mat          <- matrix(rnorm(G*nrow(data)),nrow(data),G)
      # not very satisfying that at least one of the component is zero...
      mat          <- mat-apply(mat,1,min)
      mat          <- mat/apply(mat,1,sum)
   }
   else
      mat          <- matrix(1,nrow(data),G)
   #
   msEst        <- mstep(modelName = modelName, data=data,z=mat)
   # an initial model is then estimated and provided as starting point to the em-algorithm
   em_model     <- em(modelName = msEst$modelName, data = data,parameters = msEst$parameters)
   return(em_model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
