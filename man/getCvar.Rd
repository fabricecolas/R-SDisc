\name{getCvar}
\alias{getCvar}
\title{ Returns the variables involved in cluster analysis}
\description{ Given the settings of the cluster analysis, returns the different
variable on which we should confine our cluster analysis.}
\usage{ getCvar(x) }
\arguments{
  \item{x}{ A character matrix specifying the settings of a cluster analysis, see \code{setCdata}.}
}
\value{ Returns a character vector naming the different variables.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{setCdata}}, \code{\link{genCdataSettings}},}
