\name{compareCresult}
\alias{compareCresult}
\title{ Compare the cluster results of one analysis}
\description{ This function takes one or two cluster analysis (\code{cresult})
as parameter and compare them. If a single analysis is provided as argument,
then the top models are cross-compared. Else, if two different analysis are
provided, then the query parameter is used to tell which model to
cross-compare. 

The function returns a list of cross-comparison tables. But it also saves into
CSV files its output.  }
\usage{ compareCresult(x, y, probs = seq(0, 1, 0.1)) }
\arguments{
  \item{x}{ a cresult,}
  \item{y}{ a cresult,}
  \item{probs}{ quantiles.}
}
\author{ Fabrice Colas }
