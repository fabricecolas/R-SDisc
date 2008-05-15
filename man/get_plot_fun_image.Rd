\name{get_plot_fun_image}
\alias{get_plot_fun_image}
\title{ Unexecuted plotting function for an image}
\description{ This function returns a function that can plot an image upon its
execution.  The matrix, i.e. usually the center of each cluster, is passed as
the only necessary plot argument whereas the others are pre-set. The pre-sets
are those of \code{canalysis_config}.  }
\usage{ get_plot_fun_image(def_visu_params, def_title_prefix) }
\arguments{
  \item{def_visu_params}{ Default value, taken from \code{attr(canalysis_config,"visu_params")}  }
  \item{def_title_prefix}{ Default value, taken from \code{attr(canalysis_config,"prefix")}   }  
}
\value{ Returns a function able to plot an image when called.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{get_plot_fun_dendro}}, \code{\link{get_plot_fun_image}},
\code{\link{get_plot_fun_legend}}, \code{\link{get_plot_fun_parcoord}}}
% \keyword{ ~kwd1 }
