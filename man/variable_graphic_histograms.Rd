\name{variable_graphic_histograms}
\alias{variable_graphic_histograms}
\title{Plot the histograms for each variables of a data set}
\description{ Function that returns histogram plots for each variables of a
data set. When available in the \code{canalysis_data} structure, it eventually
overload the histograms with additional information on the transformation such
as the variable empirical mean, its standard deviation or the formula for the
adjustment.  }
\usage{ variable_graphic_histograms(data, which_data = NULL) }
\arguments{
  \item{data}{ A data set of type \code{canalysis_data},}
  \item{which_data}{ If null, take the transformed data, otherwise take the
  original data of \code{canalysis_data}}
}
\details{ The plotting outputs its graphics on the current device.  }
\references{ }
\author{ Fabrice Colas }
\seealso{ \code{\link{variable_graphic_report}},\code{\link{variable_graphic_histograms}}}
\examples{ }
% __ONLY ONE__ keyword per line
