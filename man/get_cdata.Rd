\name{get_cdata}
\alias{get_cdata}
\title{ Returns the data matrix for the cluster analysis}
\description{ Given a data set appropriately transformed of type
\code{cdata}, this function returns the data matrix for use in a
cluster analysis.  }
\usage{ get_cdata(data, which_data = NULL) }
\arguments{
  \item{data}{ A data set of type \code{cdata} as prepared by \code{set_cdata},}
  \item{which_data}{ A character element with values NULL or \code{data_o}.}
}
\value{ A data matrix that we restrict to the set of variables involved in the
cluster analysis.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{get_canalysis_variables}}, \code{\link{set_cdata}}}
\keyword{}
