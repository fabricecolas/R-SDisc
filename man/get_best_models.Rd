\name{get_best_models}
\alias{get_best_models}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get the list of top ranking models in a cresult }
\description{ Returns a character vector of the type "VVV,3" that lists the top
ranking models in the \code{cresult}. }
\usage{ get_best_models(cresult) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cresult}{ A cluster analysis data structure }
}
\value{ A character vector with elements of the type "modelName, number of
cluster".  }
\author{ Fabrice Colas }
\keyword{}
