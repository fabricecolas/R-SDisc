\name{compare_cresult}
\alias{compare_cresult}
\title{ Compare the cluster results of one analysis}
\description{ This function takes one or two cluster analysis (\code{cresult})
as parameter and compare them. If a single
analysis is provided as argument, then the top models are cross-compared. Else,
if two different analysis are provided, then the query parameter is used to
tell which model to cross-compare. 

The function returns a list of cross-comparison tables. But it also saves into
CSV files its output.  }
\usage{ compare_cresult(x, query = NULL) }
\arguments{
  \item{x}{ A cluster analysis result as provided by the function \code{analysis},}
  \item{query }{ A query of the type (model name, number of clusters) with
  model name taking values in EII, VII, EEI, VEI, EVI, VVI, EEE, EEV, VEV,
  VVV.}
}
\value{
  \item{table }{ A matrix where two models are cross-compared}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{cross_compare_models}} }
