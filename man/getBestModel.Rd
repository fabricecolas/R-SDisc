\name{getBestModel}
\alias{getBestModel}
\title{ Get the list of top ranking models in a cresult }
\description{ Returns a character vector of the type "VVV,3" that lists the top
ranking models in the \code{cresult}. }
\usage{ getBestModel(x, n = NULL) }
\arguments{
  \item{x}{ A cluster analysis data structure of type (\code{cresult}),} 
  \item{n}{ Last rank to return.}
  }
\value{ A character vector with elements of the type "modelName, number of
cluster".  }
\author{ Fabrice Colas }
