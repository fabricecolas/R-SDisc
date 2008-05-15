\name{compare_canalysis}
\alias{compare_canalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compare the cluster results of one analysis}
\description{
This function takes one or more cluster analysis as parameter and it conducts a
comparative analysis of its models. If a single analysis is provided as
argument, then the top models are cross-compared. Else, if two different
analysis are provided, then the query parameter is used to tell which model to
cross-compare. 

The function returns a list of cross-comparison tables. But it also saves into
CSV files its output.
}
\usage{
compare_canalysis(analysis1, analysis2 = NULL, query = NULL, config)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{analysis1 }{ A cluster analysis result as provided by the function 'analysis',}
  \item{analysis2 }{ Either NULL or a cluster analysis result as provided by the function 'analysis',}
  \item{query }{ A query of the type (model name, number of clusters) with model name taking values in EII, VII, EEI, VEI, EVI, VVI, EEE, EEV, VEV, VVV,}
  \item{config }{ A general configuration of type 'canalysis_config'. }
}
\details{ }
\value{
  \item{table }{ A matrix where two models are cross-compared}
}
\references{ }
\author{ Fabrice Colas }
\note{  ~Make other sections like Warning with \section{Warning }{....} ~ }
\seealso{ \code{\link{cross_compare_models}} }
\examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
