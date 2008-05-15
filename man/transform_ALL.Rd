\name{transform_ALL}
\alias{transform_ALL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Routine function that transforms data given its arguments}
\description{ Routine function called by 
\code{transform_ABSMAX}, \code{transform_AVG}, \code{transform_L1},
\code{transform_MAX}, \code{transform_MEDIAN}, \code{transform_MIN} and
\code{transform_SIGMA} to transform the data provided provided two operations.
The first operation computes the estimate given the chosen transformation and
the second applies the estimate to the data. For instance, when normalizing,
the apply function is a division and when centering around the mean or the
median, the operation is the minus.  

The function outputs a transformed vector but also, for reproducibility on new
data, a \code{tdata} structure which stores the transformation properties. So,
when a non-null \code{tdata} structure is provided, the function transforms
\code{var} with respect to the value defined in \code{tdata}.  }
\usage{ transform_ALL(data 
        , tfun = list(  function(x) { return(0) }, 
                        function(x){ return(do.call("+", x)) })
        , tdata) }
\arguments{
  \item{data}{ A data structure of type \code{canalysis_data} as constructed by
               \code{set_canalysis_data},}
  \item{tfun}{ A list of two functions. One returns its transformation-estimate
               given a numeric vector. And the other applies the transformation
               to a numeric vector,}
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
%\keyword{ ~kwd1 }
