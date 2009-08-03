\name{plot.cdata}
\alias{plot.cdata}
\title{ Generic function characterizing graphically a data set}
\description{ For objects of type (\code{cdata}), this generic
plotting function reports boxplots and histograms for each variable of the data
set. On the histograms, additional textual information logs the result of the
transformation for each variable. For instance, the estimated mean or median or
the standard deviation might be reported.}
\usage{ \method{plot}{cdata}(x, ...) }
\arguments{
   \item{x}{ An object of type \code{cdata},}
   \item{...}{...}
}
\author{ Fabrice Colas }
\seealso{\code{\link{getPlotFun}}, \code{\link{setCdata}}}
