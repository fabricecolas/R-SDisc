\name{write.cresult}
\alias{write.cresult}
\title{ Generic print function that writes cluster result models as CSV}
\description{ For objects of type (\code{cresult}), this 
function writes model files as CSV. Possibly, a query parameter select a
particular model to write but by default, the top ranking ones are selected.
The output models contain the cluster membership probabilities for each
observation as well as the most likely cluster affectation.}
\usage{ write.cresult(x, query=NULL) }
\arguments{
   \item{x}{ An object of type \code{cresult},}
   \item{query}{ A valid reference (\code{query}) to the particular cluster
      result to illustrate. A reference can be of two types: an integer in
      \code{1:length(cresult)}, or a character string taken in
      \code{names(cresult)}. }
}
\author{ Fabrice Colas }
\seealso{\code{\link{cross_compare_models}}, \code{\link{set_cresult}}}
\keyword{}
