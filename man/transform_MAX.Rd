\name{transform_MAX}
\alias{transform_MAX}
\title{Normalize a variable relative to its L2 norm}
\description{ This function normalizes the values of variable \code{var} in
\code{data} relatively to the maximum value. 

The function outputs a transformed vector but also, for reproducibility on new
data, a \code{tdata} structure which stores the transformation properties. So,
when a non-null \code{tdata} structure is provided, the function normalizes
\code{var} with respect to the value defined in \code{tdata}.  }
\usage{ transform_MAX(data, tdata) }
\arguments{
  \item{data}{ A data structure of type \code{canalysis_data} as constructed by
               \code{set_canalysis_data},}
  \item{tdata}{ A transformation-data structure of type \code{tdata} as
               constructed by \code{set_tdata}.}
}
\details{ }
\value{
  It returns a list made of two elements:
  \item{vdata }{ A numeric vector of length \code{nrow(attr(data,"data_o"))}, eventually exhibiting NA's,}
  \item{tdata }{ A structure of type \code{tdata} where the transformation properties are defined.}
}
\references{ }
\author{ Fabrice Colas }
\note{ }
\seealso{ \code{\link{transform_ABSMAX}}, \code{\link{transform_ALL}}, \code{\link{transform_AVG}}, \code{\link{transform_L1}}, \code{\link{transform_MAX}}, \code{\link{transform_MEDIAN}}, \code{\link{transform_MIN}}, \code{\link{transform_SIGMA}}, \code{\link{transform_adjust}}, \code{\link{transform_canalysis_data}}, \code{\link{transform_sibordering}}}
\examples{ }
%\keyword{ ~kwd1 }%- Also NEED an '\alias' for EACH other topic documented here.
