\name{set_canalysis_config}
\alias{set_canalysis_config}
\title{ Prepare a cluster analysis configuration object }
\description{ This function contructs an object defining the global
settings of the cluster analysis to run.  }
\usage{ set_canalysis_config(prefix = NULL, visu_params = NULL, stats_fun =
list(oddratios = get_stats_fun(fun_name = "oddratios")), bbox_threshold = 5) }
\arguments{
  \item{prefix}{ Prefix when generating file outputs, automatically add the date of the day in front,}
  \item{visu_params}{ List containing visualization parameters used by the
  different graphic function, namely \code{legend_on} (name of the sub plot
  having the legend. It should match a regular group in
  \code{attr(canalysis_data,"settings")[,"visu_groups"]}),
  \code{legend_xcoord} (the x-axis coordinate of the legend),
  \code{legend_ycoord} (the y-axis coordinate of the legend), \code{cex} (the
  regular character expansion parameter), \code{color_gradient} (an hexadecimal
  color gradient in use for the heatmaps), \code{xlim=list(min,max)} and
  \code{ylim=list(min,max)} (the x and y axis boundaries),}
  \item{stats_fun}{ A list of functions that will be executed to derive additional
  statistics on the cluster results. These function calls must be accessed
  through \code{get_stats_fun()} that accepts as parameter the name of the
  function \code{fun_name} (oddratios, gen_naive_bayes, gen_knn, gen_svm,
  auuc) and eventually other settings arguments.  }
  \item{bbox_threshold}{ ~~Describe \code{bbox_threshold} here~~ }
}
\value{
  Return a data structure: 
  \item{visu_params }{Visualization parameters,}
  \item{prefix }{Generated prefix for the current analysis,}
  \item{stats_fun }{The list of functions to compute the statistics}
  \item{bbox_threshold }{The threshold for the BIC score, e.g. 5\%}
  \item{class }{'canalysis_config'}
}
\references{ }
\author{ Fabrice Colas }
\seealso{ \code{\link{set_calgo}}, \code{\link{set_canalysis_config}},
\code{\link{set_canalysis_data}}, \code{\link{set_canalysis_result}},
\code{\link{set_plot_data}}, \code{\link{set_tdata}} }
%\keyword{ ~kwd1 }
