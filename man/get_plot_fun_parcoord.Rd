\name{get_plot_fun_parcoord}
\alias{get_plot_fun_parcoord}
\title{ Unexecuted plotting function for a parallel coordinate sub plot}
\description{ This function returns a function that can do a parallel
coordinates sub plot upon its execution.  The pattern to plot is passed as its
only necessary argument whereas the others are pre-set. The pre-sets are those
of \code{canalysis_config}.  }
\usage{ get_plot_fun_parcoord(def_visu_params, def_color_sel, def_T_s,
def_title_prefix) }
\arguments{
  \item{def_visu_params}{ Default value, taken from \code{attr(canalysis_config,"visu_params")}  }
  \item{def_color_sel}{Default value, taken from \code{get_coloring_scheme()}}
  \item{def_T_s}{ Default value, taken from \code{attr(canalysis_data,"settings")[,"visu_ycoord"]}}
  \item{def_title_prefix}{ Default value, taken from \code{attr(canalysis_config,"prefix")}   }  
}
\value{ Returns a function able to do a parallel coordinates plot when called.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{get_plot_fun_dendro}}, \code{\link{get_plot_fun_image}},
\code{\link{get_plot_fun_legend}}, \code{\link{get_plot_fun_parcoord}}}
% \keyword{ ~kwd1 }
