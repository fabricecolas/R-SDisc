\name{mclust2}
\alias{mclust2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mclust2(data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = NULL, warn = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{G}{ ~~Describe \code{G} here~~ }
  \item{modelNames}{ ~~Describe \code{modelNames} here~~ }
  \item{prior}{ ~~Describe \code{prior} here~~ }
  \item{control}{ ~~Describe \code{control} here~~ }
  \item{initialization}{ ~~Describe \code{initialization} here~~ }
  \item{warn}{ ~~Describe \code{warn} here~~ }
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
function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = NULL, warn = FALSE, ...) 
{
    mc                  <- match.call(expand.dots = FALSE)
    mc[[1]]             <- as.name("mclust2_BIC")
    Bic                 <- eval(mc, parent.frame())
    G                   <- attr(Bic, "G")
    modelNames          <- attr(Bic, "modelNames")
    Sumry               <- summary_mclust2_BIC(Bic, data, G = G, modelNames = modelNames)
    if (!(length(G) == 1)) {
        bestG           <- length(unique(Sumry$cl))
        if (bestG == max(G)) 
            warning("optimal number of clusters occurs at max choice")
        else if (bestG == min(G)) 
            warning("optimal number of clusters occurs at min choice")
    }
    attr(Bic, "n")      <- attr(Bic, "warn") <- NULL
    attr(Bic, "initialization") <- attr(Bic, "control") <- NULL
    attr(Bic, "d")      <- attr(Bic, "returnCodes") <- attr(Bic, "class") <- NULL
    oldClass(Sumry)     <- NULL
    Sumry$bic           <- Sumry$bic[1]
    ans                 <- c(list(out = Bic), Sumry)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty","out")
    structure(ans[orderedNames], class = "mclust2_mclust2_Mclust")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
