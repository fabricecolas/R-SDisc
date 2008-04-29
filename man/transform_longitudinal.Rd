\name{transform_longitudinal}
\alias{transform_longitudinal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_longitudinal(data, fun, canalysis_variables = parkinson_canalysis_variables, effect = "duration")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{fun}{ ~~Describe \code{fun} here~~ }
  \item{canalysis_variables}{ ~~Describe \code{canalysis_variables} here~~ }
  \item{effect}{ ~~Describe \code{effect} here~~ }
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
function(data, fun, canalysis_variables = parkinson_canalysis_variables,effect="duration"){
   data_out             <- data
   m                    <- fun(as.data.frame(data[,c(canalysis_variables,effect),1]))
   m_effect             <- 
   data_out[row.names(m[["data"]]),colnames(m[["data"]]),1]     <- data.matrix(m[["data"]])
   for(i in 2:dim(data)[[3]]){
      # WE WANT TO REMOVE THE EFFECT OF DISEASE DURATION (AKA THE TIME) FOR BASELINE,
      # BUT WE PRESERVE THE TIME EFFECT IN THE LONGITUDINAL ANALYSIS. 
      data_tmp          <- as.data.frame(cbind(data[,canalysis_variables,i],data[,effect,1]))
      colnames(data_tmp)[ncol(data_tmp)] <- effect 
      m_tmp          <- fun(data_tmp,model = m[["model"]])
      data_out[row.names(m_tmp[["data"]]),canalysis_variables,i]  <- data.matrix(m_tmp[["data"]])[,canalysis_variables]
      }
   return(data_out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
