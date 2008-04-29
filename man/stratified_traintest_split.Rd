\name{stratified_traintest_split}
\alias{stratified_traintest_split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stratified_traintest_split(data, K = 10, perc = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{perc}{ ~~Describe \code{perc} here~~ }
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
function(data,K=10,perc=0.7){
   indexes <- list()
   G <- length(table(data$class))
   for(g in 1:G){
      strata_idx <- nrow(data[data$class == g,])
      indexes[[g]] <- matrix(0,0,strata_idx * (1-perc)) 
      set.seed(6013)
      for(k in 1:K)
            indexes[[g]] <- rbind(indexes[[g]],sample(1:strata_idx)[1:(strata_idx * (1-perc))])
      }
   return(list(indexes=indexes,K=K))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
