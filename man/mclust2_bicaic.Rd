\name{mclust2_bicaic}
\alias{mclust2_bicaic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mclust2_bicaic(modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modelName}{ ~~Describe \code{modelName} here~~ }
  \item{loglik}{ ~~Describe \code{loglik} here~~ }
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{d}{ ~~Describe \code{d} here~~ }
  \item{G}{ ~~Describe \code{G} here~~ }
  \item{noise}{ ~~Describe \code{noise} here~~ }
  \item{equalPro}{ ~~Describe \code{equalPro} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function (modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE,...)
{
   modelName    <- switch(EXPR = modelName, XII = "EII", XXI = "EEI",XXX = "EEE", modelName)
   if (G == 0) {
      if (!noise)
      stop("undefined model")
      nparams <- 1
   }
   else {
         nparams <- nVarParams(modelName, d, G) + G * d
      if (!equalPro)
         nparams <- nparams + (G - 1)
      if (noise)
         nparams <- nparams + 2
   }
   #
   # RETURN VECTOR WITH LOG BIC AND AIC SCORES
   #
   return(c(2 * loglik - nparams * logb(n), 2 * loglik - 2 * nparams))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
