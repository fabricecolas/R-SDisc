\name{funStatsOddratios}
\alias{funStatsOddratios}
\title{ For a factor, odd of the cluster as compared to the data set distribution}
\description{ For each group of variables, i.e. factors, we compute the odd of
a cluster distribution as compared to the one of the dataset.  }
\usage{ funStatsOddratios(data, class, fun_midthreshold = median) }
\arguments{
  \item{data}{ A data matrix,}
  \item{class}{ A class label,}
  \item{fun_midthreshold}{ A function (\code{mean}, \code{median}) to compute
  the threshold.}
}
\author{ Fabrice Colas }
