\name{em_clustering}
\alias{em_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
em_clustering(data, G_set = 1:9, modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV"), rseed = 6013)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{modelNames}{ ~~Describe \code{modelNames} here~~ }
  \item{rseed}{ ~~Describe \code{rseed} here~~ }
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
function(data, G_set=1:9,
        modelNames=c("EII", "VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"
        ),rseed=6013) {
   set.seed(rseed)
   # clustering results
   cresult_tmp      <- list()
   idx          <- 1
   BIC          <- matrix(NA,length(G_set),length(modelNames),dimnames=list(G_set, modelNames))
   for(g in G_set){
      for(m in modelNames){
         cresult_tmp[[idx]]             <- em_mbc(data=data,modelName=m, G=g)
         BIC[as.character(g),m]         <- bic(modelName = m, loglik = cresult_tmp[[idx]]$loglik,n=nrow(data),d=ncol(data),G=g)
         idx                            <- idx+1
      }
   }
   best_model_idx       <- which(BIC == max(BIC,na.rm=TRUE))[1]
   cresult              <- cresult_tmp[[best_model_idx]]
   cresult[["out"]]     <- cresult_tmp
   attr(cresult[["out"]],"BIC") <- BIC
   return(cresult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
