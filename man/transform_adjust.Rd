\name{transform_adjust}
\alias{transform_adjust}
\title{Adjust a variable relative for a known effect like the time}
\description{ Reduces the variability of the variable \code{var} to the one
which is not explained by the time-variable for example. The formula provided
as argument is passed to the \code{lm} function and therefore, any valid
formula on the data is accepted. 

The function outputs a transformed vector but also, for reproducibility on new
data, a \code{tdata} structure which stores the transformation properties. So,
when a non-null \code{tdata} structure is provided, the function adjust 
\code{var} with respect to the value defined in \code{tdata}.  }
\usage{ transform_adjust(data, tdata, tformula) }
\arguments{
  \item{data}{ A data structure of type \code{cdata} as constructed by
               \code{set_cdata},}
  \item{tdata}{ A transformation-data structure of type \code{tdata} as
               constructed by \code{set_tdata}.}
  \item{tformula}{ A valid \code{lm} formula defining the adjustment.}
}
\value{
  It returns a list made of two elements:
  \item{vdata }{ A numeric vector of length \code{nrow(attr(data,"data_o"))}, eventually exhibiting NA's,}
  \item{tdata }{ A structure of type \code{tdata} where the transformation properties are defined.}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{transform_adjust}}, \code{\link{transform_ABSMAX}},
\code{\link{transform_ALL}}, \code{\link{transform_AVG}},
\code{\link{transform_L1}}, \code{\link{transform_MAX}},
\code{\link{transform_MEDIAN}}, \code{\link{transform_MIN}},
\code{\link{transform_SIGMA}}, \code{\link{transform_adjust}},
\code{\link{transform_cdata}}, \code{\link{transform_sibordering}}}
\keyword{}
