\name{cresult_graphic_characteristics}
\alias{cresult_graphic_characteristics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate graphic summaries that characterize the cluster results}
\description{ 
This function produces graphics to characterize the main properties of each
cluster. It uses parallel coordinates visualization paradigm and colored
image of the cluster average patterns, dendrograms on the clusters and on
the variables.}
\usage{ cresult_graphic_characteristics(cresult, canalysis_data,
canalysis_config) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cresult}{ A cluster result as provided by \code{analysis} or \code{calgo},}
  \item{canalysis_data}{ A data container for cluster analysis of type
      \code{canalysis_data}, as constructed by \code{set_canalysis_data}}
  \item{canalysis_config}{ A cluster analysis general configuration of type
      \code{canalysis_config}, as constructed by \code{set_canalysis_config}}
}
\details{ This function prepares and execute the R code for graphic display.
The prepared graphic function are also stored back in the cluster result
(\code{cresult}) to eventually plot again parts of the cluster analysis. }
\value{
  \item{cresult }{ A cluster result as provided by \code{analysis} or \code{calgo}.}
}
\references{ }
\author{ Fabrice Colas }
\note{ ~Make other sections like Warning with \section{Warning }{....} ~ }
\seealso{ \code{\link{cresult_numeric_characteristics}}, \code{\link{compare_canalysis}}, \code{\link{analysis}}}
\examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
