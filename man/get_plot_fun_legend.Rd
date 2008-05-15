\name{get_plot_fun_legend}
\alias{get_plot_fun_legend}
\title{ Unexecuted plotting function for a legend}
\description{ This function returns a function that can plot a legend upon its
execution.  The ordered class counts are passed as its only necessary argument
whereas the others are pre-set. The pre-sets are those of
\code{canalysis_config}.  }
\usage{ get_plot_fun_legend(def_visu_params, def_color_settings) }
\arguments{
  \item{def_visu_params}{ Default value, taken from \code{attr(canalysis_config,"visu_params")}  }
  \item{def_color_settings}{Default value, taken from \code{get_coloring_scheme()}}
}
\value{ Returns a function able to plot a legend when called.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{get_coloring_scheme}}, \code{\link{get_plot_fun_dendro}}, \code{\link{get_plot_fun_image}},
\code{\link{get_plot_fun_legend}}, \code{\link{get_plot_fun_parcoord}}}
% \keyword{ ~kwd1 }
