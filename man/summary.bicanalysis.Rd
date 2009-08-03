\name{summary.bicanalysis}
\alias{summary.bicanalysis}
\title{ Generic function to summarize a bicanalysis}
\description{ For objects of type (\code{bicanalysis}), this generic
summary function reports a number of BIC tables for model selection from a
subtype discovery analysis (mean, median, ranking and confidence interval). }
\usage{ \method{summary}{bicanalysis}(object, fun = "mean", bic =
"relativeBic", type = "latex", caption = NULL, label = "tab:bic", fmt = "\%.2f",
...)
}
\arguments{
   \item{object}{ An object of type \code{bicanalysis},}
   \item{fun}{ in \code{mean}, \code{median} or \code{rank},}
   \item{bic}{ in \code{relativeBic} or \code{BIC},}
   \item{type}{ in \code{html}, \code{latex} or else, which by default will be \code{term},}
   \item{caption}{ a caption in the case of an \code{html} or \code{latex} output, the
      default behaviour is to concatenate \code{fun} and \code{bic},}
   \item{label}{ The label for cross-reference in latex documents, }
   \item{fmt}{ a type of number formatting (see \code{fmt} parameter of
      \code{?sprintf}),}
   \item{...}{...}
}
\author{ Fabrice Colas }
