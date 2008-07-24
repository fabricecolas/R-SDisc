\name{get_canalysis_variables}
\alias{get_canalysis_variables}
\title{ Returns the variables involved in cluster analysis}
\description{ Given the settings of the cluster analysis, returns the different
variable on which we should confine our cluster analysis.}
\usage{ get_canalysis_variables(settings) }
\arguments{
  \item{settings}{ A character matrix specifying the settings of a cluster analysis, see \code{set_cdata}.}
}
\value{ Returns a character vector naming the different variables.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{set_cdata}}, \code{\link{generate_cdata_settings}},}
\keyword{}
