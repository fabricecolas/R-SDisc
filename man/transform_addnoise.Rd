\name{transform_addnoise}
\alias{transform_addnoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_addnoise(data, model = NULL, canalysis_variables = NULL, type = "addnoise", relative_noise_amount = 0.01, rand_seed = 6013)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{canalysis_variables}{ ~~Describe \code{canalysis_variables} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{relative_noise_amount}{ ~~Describe \code{relative_noise_amount} here~~ }
  \item{rand_seed}{ ~~Describe \code{rand_seed} here~~ }
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
function(
        data,
        model=NULL,
        canalysis_variables = NULL,
        type="addnoise",
        relative_noise_amount = 0.01,
        rand_seed = 6013
        )
{
   if(is.null(model)){
      for(o in which(canalysis_variables \%in\% colnames(data))){
         local_name             <- colnames(data)[o]
         local_rand_seed        <- rand_seed+o
         set.seed(local_rand_seed)
         local_rand_vector      <- data.frame(rnorm(    length(data[,o])
                                                        , mean = mean(data[,o])
                                                        , sd = relative_noise_amount * sd(data[,o])
                                                   )
                                             , row.names = row.names(data))
         colnames(local_rand_vector) <- "values"
         model[[local_name]]    <- list(  rand_seed = local_rand_seed
                                        , rand_vector = local_rand_vector
                                        , relative_noise_amount = relative_noise_amount)
         # (EXPECT NO NA'S WHILE SUMMING BECAUSE "0.54345<-NA+0.54345"...)
         data[,o] <- data[,o] + local_rand_vector[,1]
      }
   }
   else{
      for(o in colnames(data))
         data[,o] <- data[,o] + model[[o]][["rand_vector"]][row.names(data),]
   }
   return(list(data=data,type=type,model=model))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
