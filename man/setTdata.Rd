\name{setTdata}
\alias{setTdata}
\title{ Prepare a transform data object}
\description{ This function returns a data structure (\code{tdata}) which
stores the parameters of a particular transform. So, given new data, we use
this information to transform them according to previously estimated transform
parameters.  }
\usage{ setTdata(testimate = NULL, var = NULL) }
\arguments{
  \item{testimate}{ The transformation estimate(s),}
  \item{var}{ The variable which should be transformed. It is especially used
  with \code{transform_adjust}.}
}
\value{ Returns a \code{tdata} structure.}
\author{ Fabrice Colas }
\seealso{ \code{\link{setCdata}}, \code{\link{setCresult}}, }
