\name{transform_variables}
\alias{transform_variables}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_variables(data, fun_transform = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{fun_transform}{ ~~Describe \code{fun_transform} here~~ }
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
function(data,fun_transform = NULL){
   #
   # APPLY USER-PROVIDED DATA TRANSFORMATION
   #
   # data_transformed <- data[["cc"]][,unlist(data[["sumscore_groups"]])]
   data_transformed <- data[["cc"]]
   model_out <- list()
   for(fun_name in names(fun_transform)){
      model_out[[fun_name]] <- fun_transform[[fun_name]](data_transformed)
      data_transformed      <- model_out[[fun_name]][["data"]]
   }
   return(list(data=data_transformed,model=model_out))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
