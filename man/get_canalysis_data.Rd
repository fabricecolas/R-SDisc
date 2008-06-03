\name{get_canalysis_data}
\alias{get_canalysis_data}
\title{ Returns the data matrix for the cluster analysis}
\description{ Given a data set appropriately transformed of type
\code{canalysis_data}, this function returns the data matrix for use in a
cluster analysis.  }
\usage{ get_canalysis_data(data, which_data = NULL) }
\arguments{
  \item{data}{ A data set of type \code{canalysis_data} as prepared by \code{set_canalysis_data},}
  \item{which_data}{ A character element with values NULL or \code{data_o}.}
}
\value{ A data matrix that we restrict to the set of variables involved in the
cluster analysis.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{get_canalysis_variables}}, \code{\link{set_canalysis_data}}}
\keyword{}
