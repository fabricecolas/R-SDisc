\name{print.cresult}
\alias{print.cresult}
\title{ Generic print function with the characteristics of a cluster result}
\description{ For objects of type (\code{cresult}), this generic
print function reports cross comparison of the most likely cluster results
(\code{compareModels}).}
\usage{ \method{print}{cresult}(x, m1 = NULL, m2 = NULL, type = "term", label = "tab:jointdistrib",...) }
\arguments{
   \item{x}{ An object of type \code{cresult},}
   \item{m1}{ A model name as returned by \code{getBestModel},}
   \item{m2}{ Same than \code{m2},}
   \item{type}{ The type of output to generate from XTable (html, latex, or term),}
   \item{label}{ The label for cross-referencing in latex,}
   \item{...}{...}
}
\author{ Fabrice Colas }
