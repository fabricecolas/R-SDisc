\name{print.bicanalysis}
\alias{print.bicanalysis}
\title{Returns as a list a BIC analysis on a set of cluster results}
\description{ This functions turns into a list a \code{bicanalysis} data
structure which was computed on a set of cluster results (\code{cresult}) }
\usage{ \method{print}{bicanalysis}(x, silent = TRUE, ...) }
\arguments{
  \item{x}{ A \code{bicanalysis} data structure, normally stored into
      \code{attr(*,'bicanalysis')} of a \code{cresult},}
  \item{silent}{ Whether we should output something, by default not,}
  \item{\dots}{ Additional arguments}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{set_bicanalysis}}, \code{\link{set_cresult}} }
