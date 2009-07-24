\name{get_plot_fun2}
\alias{get_plot_fun2}
\title{ Returns a generic unevaluated plotting function}
\description{ This function returns a function that plots a legend, a
dendrogram, or an image. The main plot settings are pre-set at the
initialization and at evaluation, it takes a single \code{cmodel} to plot the
patterns.}
\usage{ get_plot_fun2(type, cdata, title, xlim, xy, cex, pattern,
color_gradient) }
\arguments{
   \item{type}{ The type of plot function with values in (\code{plot_parcoord},
      \code{plot_image}, \code{plot_legend}, \code{plot_dendro_var},
      \code{plot_dendro_cluster}),}
   \item{cdata}{ A data set (\code{set_cdata}) to configure the plotting
      function of type \code{plot_parcoord} where several sub plots
      (\code{visu_groups}) are generated,}
   \item{title}{ A character string defining the title for the sub plot,}
   \item{xlim}{ The x-axis limits in a plot (\code{plot_parcoord}) with default values set to \code{c(-3,3)},}
   \item{xy}{ The x- and y-axis coordinates of the lower-left legend element, with default value \code{c(-2.2,0)},}
   \item{cex}{ A character expansion factor with default value set to \code{0.9},}
   \item{pattern}{ A character vector defining the type of pattern to report in
   the plots. The character vector may have a single element (\code{plot_image},
      \code{plot_parcoord}) or more (\code{plot_parcoord}). The set of possible
      patterns refers to those defined by \code{fun_pattern} in
      \code{set_cresult}, which defaults to \code{mean, median, lowquant,
      upquant},}
   \item{color_gradient}{ A gradient of color, especially for \code{plot_image}.}
}
\value{ Returns a function taking a single \code{cmodel} argument.  }
\author{ Fabrice Colas }
\seealso{\code{\link{plot_parcoord}}, \code{\link{set_cresult}}}
