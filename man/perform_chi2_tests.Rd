\name{perform_chi2_tests}
\alias{perform_chi2_tests}
\title{Perform chi2 tests of association between each feature and a class column}
\description{Iteratively, the features are tested for association with the
class column. When features are not nominal, the parameter \code{tonominal} is
used to list those ordinal and continuous features that need to be binarized.
In that case, the values are tested to be above or below the median of the
feature's distribution.}
\usage{ perform_chi2_tests(m, tonominal=FALSE, simulate.p.value=TRUE, B=10000) }
\arguments{
   \item{m}{ The data matrix with a \code{class} column,}
   \item{tonominal}{ The features to binarize,}
   \item{simulate.p.value}{ A boolean telling whether p-values should be estimated by MonteCarlo simulation,}
   \item{B}{ An integer giving the number of repeat of the MonteCarlo simulation.}
}
\author{ Fabrice Colas }
\seealso{\code{\link{association_analysis}}}
