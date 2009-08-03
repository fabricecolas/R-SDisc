\name{setBicanalysis}
\alias{setBicanalysis}
\title{Function that performs a BIC analysis on a set of cluster results}
\description{ This procedure conducts a BIC score analysis given a set of
cluster results (\code{cresult}). As we repeat modelling for clusters by
different starting points chosen at random (\code{rseed}, \code{set.seed}) this
analysis computes the mean, median, standard deviation, 2.5\% and 97.5\%
quantiles, of the BIC scores for these repeats. In addition, it also computes
different ranking to determine the most likely: model, number of clusters or
random start number.  }
\usage{ setBicanalysis(x) }
\arguments{ \item{x}{ An object of type \code{cresult}.} }
\value{
  \item{attr(*,"bicanalysis")}{In an object of type \code{cresult}, }
}
\author{ Fabrice Colas }
\seealso{ \code{\link{setCresult}} }
