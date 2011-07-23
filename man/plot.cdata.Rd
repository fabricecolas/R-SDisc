\name{plot.cdata}
\alias{plot.cdata}
\title{ Generic function characterizing graphically a data set}
\description{ For objects of type (\code{cdata}), this generic
plotting function reports boxplots and histograms for each variable of the data
set. On the histograms, additional textual information logs the result of the
transformation for each variable. For instance, the estimated mean or median or
the standard deviation might be reported.}
\usage{ plot.cdata(x, device="PS", ...) }
\arguments{
   \item{x}{ An object of type \code{cdata},}
   \item{device}{ A character string telling the type of device to output the
      plot (in \code{"PS", "PNG"}). }
   \item{...}{...}
}
\author{ Fabrice Colas }
\seealso{\code{\link{get_plot_fun}}, \code{\link{set_cdata}}}
\keyword{}
