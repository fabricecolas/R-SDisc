\name{plot_parcoord}
\alias{plot_parcoord}
\title{ Unexecuted plotting function for a parallel coordinate sub plot}
\description{ This function returns a function that can do a parallel
coordinates sub plot upon its execution.  The pattern to plot is passed as its
only necessary argument whereas the others are pre-set. The pre-sets are those
of \code{canalysis_config}.  }
\usage{ plot_parcoord(xlim = c(-3,3) , cex = NULL , title = NULL , T_s = NULL ,
pattern = NULL) }
\arguments{
   \item{xlim}{Axis limit in x with default values to \code{c(-3,3)},}
   \item{cex}{A character expension factor for the test, with default values set
   to \code{0.9}}
   \item{title}{A character string for the title of the sub plot,}
   \item{T_s}{A setting matrix which defines in particular the y-axis coordinates
   and the group of variables in each subplot \code{visu_ycoord, visu_groups},}
   \item{pattern}{ A character vector defining the different patterns to plot,
   see \code{fun_pattern} in \code{set_cresult},}
}

\value{ Returns a function able to do a parallel coordinates plot when called.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{set_cresult}}, \code{\link{get_plot_fun}}}
