\name{tSigma}
\alias{tSigma}
\title{Normalizes a variable relative to its standard deviation}
\description{ This function normalizes the values of variable \code{var} in
\code{data} relatively to the standard deviation. 

The function outputs a transformed vector but also, for reproducibility on new
data, a \code{tdata} structure which stores the transformation properties. So,
when a non-null \code{tdata} structure is provided, the function normalizes
\code{var} with respect to the value defined in \code{tdata}.  }
\usage{ tSigma(data, tdata) }
\arguments{
  \item{data}{ A data structure of type \code{cdata} as constructed by
               \code{set_cdata},}
  \item{tdata}{ A transformation-data structure of type \code{tdata} as
               constructed by \code{set_tdata}.}
}
\value{
  It returns a list made of two elements:
  \item{vdata }{ A numeric vector of length \code{nrow(attr(data,"data_o"))}, eventually exhibiting NA's,}
  \item{tdata }{ A structure of type \code{tdata} where the transformation properties are defined.}
}
\author{ Fabrice Colas }
