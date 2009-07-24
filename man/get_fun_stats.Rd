\name{get_fun_stats}
\alias{get_fun_stats}
\title{ Returns a generic unevaluated stats function}
\description{ This function returns unevaluated functions that computes
statistics on a particular cluster model (\code{cmodel}). }
\usage{ get_fun_stats(fun_name = "oddratios", ...) }
\arguments{
  \item{fun_name}{ one of \code{oddratios}, \code{gen_knn}, \code{gen_svm},
  \code{gen_naive_bayes}, \code{auuc}, \code{lambdasibs}, }
  \item{\dots}{ additional parameters to be passed to the functions.}
}
\seealso{
\code{\link{stats_logodds}},
\code{\link{stats_lambdasibs}},
\code{\link{stats_auuc}},
\code{\link{stats_generalization}},
}
\author{ Fabrice Colas }
