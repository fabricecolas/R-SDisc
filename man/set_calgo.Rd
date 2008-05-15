\name{set_calgo}
\alias{set_calgo}
\title{ Prepare a generic 'cluster alorigthm' function object}
\description{ This function prepares a cluster algorithm function call. The
subsequent execution of the function needs two arguments of type
\code{canalysis_data} and \code{canalysis_config}.}
\usage{ set_calgo(calgo = "model_based_clustering", args = list(G_set = 1:9,
model_names = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV",
"VVV"))) }
\arguments{
  \item{calgo}{ A character string naming the function to call for execution,}
  \item{args}{ A list of arguments that will be passed to this
  function in addition to the data.}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{set_calgo}}, \code{\link{set_canalysis_config}},
\code{\link{set_canalysis_data}}, \code{\link{set_canalysis_result}},
\code{\link{set_plot_data}}, \code{\link{set_tdata}} }
%\keyword{ ~kwd1 }
