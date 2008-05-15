\name{area_under_uncertainty_curve}
\alias{area_under_uncertainty_curve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Area under the quantile-curve of the clustering uncertainty}
\description{
Each observation belong to a cluster with a given certainty. And we can refer
as the clustering uncertainty the reciprocal, i.e.  1 minus the cluster
likelihood membership. Then, the function (area_under_uncertainty_curve())
returns the mean value of the clustering uncertainties. 

We can understand this quantity as the average level of uncertainty to cluster
the observation for the current model.  
}
\usage{
area_under_uncertainty_curve(class)
}
\arguments{
  \item{class}{ A data matrix of the cluster membership likelihood.}
}
\details{
  
}
\value{
  \item{auuc}{ Average value of the cluster uncertainty,}
  \item{auuc_sd}{ Standard deviation of the cluster uncertainty,}
  \item{out}{ 2-by-2 matrix defining the output-print of the table.}
}
\references{ }
\author{ Fabrice Colas }
\note{ }
\seealso{  \code{\link{get_stats_fun}, \link{evaluate_generalization},\link{statistics_logodds}}}
\examples{ }
