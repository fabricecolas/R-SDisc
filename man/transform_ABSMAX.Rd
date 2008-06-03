\name{transform_ABSMAX}
\alias{transform_ABSMAX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalize a variable relatively to its absolute maximum}
\description{ This function normalizes the values of variable \code{var} in
\code{data} relatively to its absolute maximum. 

The function outputs a transformed vector but also, for reproducibility on new
data, a \code{tdata} structure which stores the transformation properties. So,
when a non-null \code{tdata} structure is provided, the function normalizes
\code{var} with respect to the value defined in \code{tdata}.  }
\usage{ transform_ABSMAX(data, tdata) }
\arguments{
  \item{data}{ A data structure of type \code{canalysis_data} as constructed by
               \code{set_canalysis_data},}
  \item{tdata}{ A transformation-data structure of type \code{tdata} as
               constructed by \code{set_tdata}.}
}
\value{
  It returns a list made of two elements:
  \item{vdata }{ A numeric vector of length \code{nrow(attr(data,"data_o"))}, eventually exhibiting NA's,}
  \item{tdata }{ A structure of type \code{tdata} where the transformation properties are defined.}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{transform_ABSMAX}}, \code{\link{transform_ALL}},
\code{\link{transform_AVG}}, \code{\link{transform_L1}},
\code{\link{transform_MAX}}, \code{\link{transform_MEDIAN}},
\code{\link{transform_MIN}}, \code{\link{transform_SIGMA}},
\code{\link{transform_adjust}}, \code{\link{transform_canalysis_data}},
\code{\link{transform_sibordering}}}
\keyword{}
