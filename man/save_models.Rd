\name{save_models}
\alias{save_models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
save_models(data, filePrefix = NULL, threshold = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{threshold}{ ~~Describe \code{threshold} here~~ }
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
function(data,filePrefix=NULL,threshold=5){
   m_set                <- get_model_set(data[["cluster_analysis"]][["bic_table_relative"]],bbox_threshold=threshold)
   for(l in 1:length(data[["cluster_analysis"]][["out"]])){
      local_canalysis   <- data[["cluster_analysis"]][["out"]][[l]]
      g_match           <- match(local_canalysis[["G"]],m_set[,"G"])
      m_match           <- match(local_canalysis[["modelName"]],m_set[,"modelName"])
      if( !is.na(g_match) & !is.na(m_match)){
	data_out        <- cbind(         ids=row.names(data[["transformed"]])
                                        , class=map(local_canalysis[["z"]])
                                        , local_canalysis[["z"]])
        colnames(data_out)[3:ncol(data_out)] <- 1:local_canalysis[["G"]]
        fileNameCSV     <- paste2(filePrefix,"_Class_",local_canalysis[["modelName"]],"-",local_canalysis[["G"]],".csv") 
        print(fileNameCSV)
        write.csv2(data_out,file=fileNameCSV)
      }
   }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
