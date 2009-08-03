\name{getBicStats}
\alias{getBicStats}
\title{From a long BIC table, returns one in wide format by eventually aggregating the results}
\description{ From a long BIC table where are repeated for each experiment the
name of the model, the number of clusters and the random start number, it
produces a wide table by eventually aggregating the results.  }
\usage{ getBicStats(x, fun="mean", vAgg="BIC", vX="G", vY="modelName", fmt="\%.2f") }
\arguments{
  \item{x}{ A long BIC table with repeated factors, }
  \item{fun}{ A regular function to aggregate the result such as \code{mean},}
  \item{vAgg}{ Which variable is to aggregate, by default "BIC" }
  \item{vX}{ Variable to be considered for the lines (\code{var_x}),}
  \item{vY}{ Variable to be considered for the lines (\code{var_y}),}
  \item{fmt}{ A character vector of format strings, each of up to 8192 bytes (see \code{sprintf}).}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{print.bicanalysis}},  \code{\link{setBicanalysis}}}
