\name{compare_canalysis}
\alias{compare_canalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compare_canalysis(analysis1, analysis2 = NULL, query = NULL, filePrefix = NULL, bbox_threshold = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{analysis1}{ ~~Describe \code{analysis1} here~~ }
  \item{analysis2}{ ~~Describe \code{analysis2} here~~ }
  \item{query}{ ~~Describe \code{query} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
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
function(analysis1,analysis2=NULL,query=NULL,filePrefix=NULL,bbox_threshold = 5){
   canalysis1 <- analysis1[["cluster_analysis"]]
   cross_comparison_results <- list()
   model      <- list()
   if(!is.null(analysis2) && !is.null(query)){
      # WE COMPARE THE RESULT OF TWO ANALYSIS FOR A GIVEN MODEL
      canalysis2 <- analysis2[["cluster_analysis"]]
      model[[1]] <- retrieve_model(canalysis1,query)
      model[[2]] <- retrieve_model(canalysis2,query)
      if(!is.null(filePrefix))
         filePrefix <- paste2(filePrefix,"_X_comparison_",gsub(",","_",query),".csv")
      if(!is.na(model[[1]][["model"]]$loglik) && !is.na(model[[2]][["model"]]$loglik))
         cross_comparison_results[[1]] <- cross_compare_models(model[[1]],model[[2]],fileNameCSV=filePrefix)
   }
   else{
      # SELECT ALL 2-BY-2 COMPARISONS BETWEEN MODELS THAT ARE WITHIN THE RELATIVE
      # 5\% BBOX TO THE BEST BIC VALUE
      vertices <- get_vertice_set(canalysis1[["bic_table_relative"]], bbox_threshold)
      if(nrow(vertices)>=1){
         # PROCEED TO 2-BY-2 COMPARISON OVER THE SET OF VERTICES
         for(v in 1:nrow(vertices)){
            # RETRIEVE CLUSTERING Z-MAPS FROM 'CANALYSIS'
            model[[1]] <- retrieve_model(canalysis1,vertices[v,1])
            model[[2]] <- retrieve_model(canalysis1,vertices[v,2])
            if(!is.null(filePrefix))
               filePrefix <- paste2(filePrefix,"_",gsub(",","_",vertices[v,1]),"-",gsub(",","_",vertices[v,2]),".csv")
            if( !is.na(model[[1]][["model"]]$loglik) && !is.na(model[[2]][["model"]]$loglik))
               cross_comparison_results[[v]] <- cross_compare_models(model[[1]],model[[2]],fileNameCSV=filePrefix)
         }
      }
   }
   return(cross_comparison_results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
