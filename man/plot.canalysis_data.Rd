\name{plot.canalysis_data}
\alias{plot.canalysis_data}
\title{ Generic function characterizing graphically a data set}
\description{ For objects of type (\code{canalysis_data}), this generic
plotting function reports boxplots and histograms for each variable of the data
set. On the histograms, additional textual information logs the result of the
transformation for each variable. For instance, the estimated mean or median or
the standard deviation might be reported.}
\usage{ plot.canalysis_data(x,doPSOutput = TRUE, ...) }
\arguments{
   \item{x}{ An object of type \code{canalysis_data},}
   \item{doPSOutput}{ A logicial value for (\code{postscript})  telling whether the
   output should be redirected to a Postscript file. }
   \item{...}{...}
}
\author{ Fabrice Colas }
\seealso{\code{\link{get_plot_fun}}, \code{\link{set_canalysis_data}}}
\keyword{}
