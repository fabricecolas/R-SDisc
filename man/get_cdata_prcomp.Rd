\name{get_cdata_prcomp}
\alias{get_cdata_prcomp}
\title{ Get the principal components explaining 95\% of the variability in a cdata object }
\description{ This method performs a principal component analysis on the set of
features defined as \code{'in_canalysis'} (\code{attr(cdata,'settings')}). By
default, those components explaining 95\% of the variability are retrieved.
Points are further projected on those orthonormal projections.  }
\usage{ get_cdata_prcomp(x, cumproportion=0.95, scale=TRUE, center=FALSE) }
\arguments{
   \item{x}{ A \code{cdata} object embeding a data set,}
   \item{cumproportion}{ Selected the components that explain, e.g. \code{cumproportion=0.95} (95\%) of the variability,}
   \item{scale}{ A boolean to indicate whether data must be scaled,}
   \item{center}{ A boolean to indicate whether data must be center.}
}
\author{ Fabrice Colas }
\seealso{\code{\link{get_cdata_princomp}}, \code{\link{get_cdata}}}
