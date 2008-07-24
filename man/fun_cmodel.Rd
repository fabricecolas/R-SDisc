\name{fun_cmodel}
\alias{fun_cmodel}
\title{ Perform a cluster analysis under different modeling settings}
\description{ Given a prepared data set (\code{cdata}), this function
performs a cluster analysis under different modeling settings. Each setting is
defined in (\code{cresult}) and iteratively, it models the data for each type
of model. It returns a set of cluster results (\code{cresult})  whose elements
are cluster models (\code{cmodel}).}
\usage{ fun_cmodel(x) }
\arguments{
  \item{x}{ A set of cluster results (\code{cresult}) as constructed by
  \code{set_cresult}. } }
\value{ Returns a \code{cresult} data structure whose cluster results are the
estimated models \code{cmodel}.  }
\references{  }
\author{ Fabrice Colas }
\seealso{ \code{\link{set_cresult}}, \code{\link{set_cdata}} }
\keyword{}
