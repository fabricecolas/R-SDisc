\name{association_analysis}
\alias{association_analysis}
\title{Perform an association analysis between a cluster result and the features}
\description{ Given the queried model (\code{q} in \code{names(x)}), perform an
association analysis between the features and the subtype membership. The
analysis proceed by repeatedly performing the chi2 test of association on a
number of bootstrap samples. }
\usage{ association_analysis(x, q=NA, B=4) }
\arguments{
   \item{x}{ A \code{cresult},}
   \item{q}{ The model on which to perform the analysis, one of \code{names(x)},}
   \item{B}{ The number of iteration of the bootstrap.}
}
\author{ Fabrice Colas }
\seealso{\code{\link{plot.association_analysis}}}
