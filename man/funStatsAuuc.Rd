\name{funStatsAuuc}
\alias{funStatsAuuc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Area under the quantile-curve of the clustering uncertainty}
\description{
Each observation belong to a cluster with a given certainty. And we can refer
as the clustering uncertainty the reciprocal, i.e.  1 minus the cluster
likelihood membership. Then, the function (\code{funStatsAuuc})
returns the mean value of the clustering uncertainties. 

We can understand this quantity as the average level of uncertainty to cluster
the observation for the current model.  
}
\usage{ funStatsAuuc(class) }
\arguments{
  \item{class}{ A data matrix of the cluster membership likelihood.}
}
\value{
  \item{auuc}{ Average value of the cluster uncertainty,}
  \item{auuc_sd}{ Standard deviation of the cluster uncertainty,}
  \item{out}{ 2-by-2 matrix defining the output-print of the table.}
}
\author{ Fabrice Colas }
