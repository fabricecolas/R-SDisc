\name{get_vertice_set}
\alias{get_vertice_set}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
get_vertice_set(relativeBic, bbox_threshold = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{relativeBic}{ ~~Describe \code{relativeBic} here~~ }
  \item{bbox_threshold}{ ~~Describe \code{bbox_threshold} here~~ }
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
function(relativeBic, bbox_threshold=5){
   relativeBic <- list( bicTable=relativeBic, 
                        # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
                        bbox = apply(relativeBic <= max(sort(relativeBic)[1:bbox_threshold]), 2,which)
                        # TO SELECT ALL MODELS BELOW A GIVEN THRESHOLD
                        # bbox = apply(abs(relativeBic) < bbox_threshold,2,which)
                        )
   relativeBic$bbox <- relativeBic$bbox[which(lapply(relativeBic$bbox,length) >0)]
   # DERIVE VERTICES TO COMPARE WITHIN THE 5\% BBOX OF THE BIC TABLE
   relativeBic$vertices <- matrix(0,0,2,dimnames=list(list(),list("m1_g1","m2_g2")))
   for(m1 in names(relativeBic$bbox))
      for(m2 in names(relativeBic$bbox))
         if(match(m2,names(relativeBic$bbox)) >= match(m1,names(relativeBic$bbox)))
            for(g1 in relativeBic$bbox[[m1]])
               for(g2 in relativeBic$bbox[[m2]])
                  if((m1 != m2 & g2 >= g1) | (m1 == m2 & g2 > g1))
                     relativeBic$vertices <- rbind(relativeBic$vertices,c(
                                                      paste(m1,row.names(relativeBic[["bicTable"]])[g1],collapse="",sep=","),
                                                      paste(m2,row.names(relativeBic[["bicTable"]])[g2],collapse="",sep=",")))
   return(relativeBic[["vertices"]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
