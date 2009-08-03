\name{getStatsFun}
\alias{getStatsFun}
\title{ Returns a generic unevaluated stats function}
\description{ This function returns unevaluated functions that computes
statistics on a particular cluster model (\code{cmodel}). }
\usage{ getStatsFun(fun_name = "oddratios", ...) }
\arguments{
  \item{fun_name}{ one of \code{oddratios}, \code{gen_knn}, \code{gen_svm},
  \code{gen_naive_bayes}, \code{auuc}, \code{lambdasibs}, }
  \item{\dots}{ additional parameters to be passed to the functions.}
}
\author{ Fabrice Colas }
