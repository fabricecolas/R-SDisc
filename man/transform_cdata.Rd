\name{transform_cdata}
\alias{transform_cdata}
\title{ Transform data for a cluster analysis}
\description{ This function applies sequencially the different transformations
as defined in settings, and this for each variable.  }
\usage{ transform_cdata(data, settings, tdata = NULL) }
\arguments{
  \item{data}{ A data set of type \code{cdata}, }
  \item{settings}{ A character matrix with row names the name of each variable
  and for column name, at least a \code{fun_transform} column. In this column
  should be defined the sequence of transformation.}
  \item{tdata }{ When reproducing results given a particular way to transform
  and to model the data, \code{tdata} argument uses previous modeling
  estimates.}
}
\value{
  \item{data }{ A data frame defining the transformed data set,}
  \item{tdata }{ A list of \code{tdata} element sequences. The first list-level
  is for the different variables of \code{data} and the second for the ordered
  sequence of transformation, which properties are memorized in \code{tdata}
  objects.}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{transform_adjust}}, \code{\link{transform_ABSMAX}},
\code{\link{transform_ALL}}, \code{\link{transform_AVG}},
\code{\link{transform_L1}}, \code{\link{transform_MAX}},
\code{\link{transform_MEDIAN}}, \code{\link{transform_MIN}},
\code{\link{transform_SIGMA}}, \code{\link{transform_adjust}},
\code{\link{transform_cdata}}, \code{\link{transform_sibordering}}}
\keyword{}
