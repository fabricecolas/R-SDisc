\name{plot.cresult}
\alias{plot.cresult}
\title{ Generic function characterizing graphically a set of cluster results}
\description{ For objects of type (\code{cresult}), this generic plotting
function draws the visual characteristicsset of a set of cluster results. This
function accepts in particular two parameters, namely \code{postscript} and
\code{query}. The first tells whether a postscript output should be generated,
whereas the second one confines the graphical characterization to a particular
queried model.}
\usage{ plot.cresult(x, device="PS", query=NULL, img_size=c(1024, 1024),...) }
\arguments{
   \item{x}{ An object of type \code{cresult},}
   \item{device}{ A character string telling the type of device to output the
   plots (in \code{"PS", "PNG"}), }
   \item{query}{ A valid reference (\code{query}) to the particular cluster
      result to illustrate. A reference can be of two types: an integer in
      \code{1:length(cresult)}, or a character string taken in
      \code{names(cresult)}, }
   \item{img_size}{ When "PNG" images are generated, default size of the plots,}
   \item{...}{ nothing }
}
\author{ Fabrice Colas }
\seealso{\code{\link{set_cresult}}, \code{\link{get_plot_fun}}}
\keyword{}
