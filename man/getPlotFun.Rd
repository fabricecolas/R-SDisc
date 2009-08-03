\name{getPlotFun}
\alias{getPlotFun}
\title{ Returns a generic unevaluated plotting function}
\description{ This function is an interface to \code{getPlotFun2}. It returns
a function that, when called with a \code{cdata} argument, returns a plot
function which will itself take a single \code{cmodel} parameter. For
additional details, see the returned ouput of \code{getPlotFun2}.}
\usage{ getPlotFun(type="plotParcoord", title=NULL, xlim=c(-3,3),
      zlim=c(-2, 2) , xy=c(-2.2,0), cex=0.7, pattern="mean",
      color_gradient=rev(brewer.pal(9, "RdBu")), range_fv=NULL) }
\arguments{
   \item{type}{ The type of plot function with values in (\code{plotParcoord},
      \code{plot_image}, \code{plot_legend}, \code{plot_dendro_var},
      \code{plot_dendro_cluster}),}
   \item{title}{ A character string defining the title for the sub plot,}
   \item{xlim}{ The x-axis limits in a plot (\code{plotParcoord}) with default values set to \code{c(-3,3)},}
   \item{zlim}{ Indicates the maxima for breaking the color gradient of the heatmap, by default \code{c(-2,2)},}
   \item{xy}{ The x- and y-axis coordinates of the lower-left legend element, with default value \code{c(-2.2,0)},}
   \item{cex}{ A character expansion factor with default value set to \code{0.9},}
   \item{pattern}{ A character vector defining the type of pattern to report in
   the plots. The character vector may have a single element (\code{plot_image},
      \code{plotParcoord}) or more (\code{plotParcoord}). The set of possible
      patterns refers to those defined by \code{fun_pattern} in
      \code{setCresult}, which defaults to \code{mean, median, lowquant,
      upquant},}
   \item{color_gradient}{ A gradient of color, especially for \code{plot_image},}
   \item{range_fv}{ A range of features to be plotted when performing analysis on series.}
}
\value{ Returns a function taking a single \code{cmodel} argument.  }
\author{ Fabrice Colas }
\seealso{\code{\link{getPlotFun2}}, \code{\link{plotParcoord}}, \code{\link{setCresult}}}
