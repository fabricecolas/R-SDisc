\name{transform_canalysis_data}
\alias{transform_canalysis_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Transform data for a cluster analysis}
\description{ This function applies sequencially the different transformations
as defined in settings, and this for each variable.  }
\usage{ transform_canalysis_data(data, settings) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data set of type \code{canalysis_data}, }
  \item{settings}{ A character matrix with row names the name of each variable
  and for column name, at least a \code{fun_transform} column. In this column
  should be defined the sequence of transformation.}
}
\details{ }
\value{
  \item{data }{ A data frame defining the transformed data set,}
  \item{tdata }{ A list of \code{tdata} element sequences. The first list-level
  is for the different variables of \code{data} and the second for the ordered
  sequenc of transformation, which properties are memorized in \code{tdata}
  objects.}
}
\references{ }
\author{ Fabrice Colas }
\note{ }
\seealso{ \code{\link{transform_adjust}}, \code{\link{transform_ABSMAX}},
\code{\link{transform_ALL}}, \code{\link{transform_AVG}},
\code{\link{transform_L1}}, \code{\link{transform_MAX}},
\code{\link{transform_MEDIAN}}, \code{\link{transform_MIN}},
\code{\link{transform_SIGMA}}, 
\code{\link{transform_adjust}}, \code{\link{transform_canalysis_data}},
\code{\link{transform_sibordering}}}
\examples{ }
% __ONLY ONE__ keyword per line
