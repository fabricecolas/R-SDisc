\name{cresult_numeric_characteristics}
\alias{cresult_numeric_characteristics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute statistics on the cluster results}
\description{
On the different models of a cluster result (\code{cresult}), this function
computes user-defined statistics. The functions are listed into the general
cluster analysis configuration of type \code{canalysis_config}.  }
\usage{ cresult_numeric_characteristics(cresult, canalysis_data,
canalysis_config) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cresult}{ A cluster result as provided by \code{analysis} or \code{calgo},}
  \item{canalysis_data}{ A data container for cluster analysis of type
      \code{canalysis_data}, as constructed by \code{set_canalysis_data}}
  \item{canalysis_config}{ A cluster analysis general configuration of type
      \code{canalysis_config}, as constructed by \code{set_canalysis_config}}
}
\details{ }
\value{
  \item{cresult }{ A cluster result as provided by \code{analysis} or \code{calgo}.}
}
\references{ }
\author{ Fabrice Colas }
\note{ }
\seealso{  \code{\link{get_stats_fun}}, \code{\link{cresult_graphic_characteristics}},
   \code{\link{compare_canalysis}}, \code{\link{analysis}} }
\examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
