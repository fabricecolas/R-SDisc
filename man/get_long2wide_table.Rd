\name{get_long2wide_table}
\alias{get_long2wide_table}
\title{From a long BIC table, returns one in wide format by eventually
aggregating the results}
\description{ From a long BIC table where are repeated for each experiment the
name of the model, the number of clusters and the random start number, it
produces a wide table by eventually aggregating the results.  }
\usage{ get_long2wide_table(x, fun_aggregate=mean, var_aggregate="BIC",
var_x="G", var_y="modelName") }
\arguments{
  \item{x}{ A long BIC table with repeated factors, }
  \item{fun_aggregate}{ A regular function to aggregate the result such as
  \code{mean},}
  \item{var_aggregate}{ Which variable is to aggregate, by default "BIC" }
  \item{var_x}{ Variable to be considered for the lines (\code{var_x}),}
  \item{var_y}{ Variable to be considered for the lines (\code{var_y}),}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{print.bicanalysis}},  \code{\link{set_bicanalysis}}}
\keyword{ }
