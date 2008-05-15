\name{variable_graphic_report}
\alias{variable_graphic_report}
\title{ Prepare a graphical report on the data set}
\description{ Prepare a graphical report on the data set that consists of a
correlation matrix, a plot of box-and-whisker and of the histograms for all
variables.  }
\usage{ variable_graphic_report(data, which_data = NULL, prefix = NULL) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data set of type \code{canalysis_data},}
  \item{which_data}{ If null, take the transformed data, otherwise take the
  original data of \code{canalysis_data}}
  \item{prefix}{ A character string for the files to output. It should be
  stored in the \code{canalysis_config}.}
}
\details{ The plotting outputs its graphics on the current device.  }
\references{ }
\author{ Fabrice Colas }
\seealso{ \code{\link{variable_graphic_report}},\code{\link{variable_graphic_histograms}}}
\examples{ }
% __ONLY ONE__ keyword per line
