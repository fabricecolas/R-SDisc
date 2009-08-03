\name{summary.cmodel}
\alias{summary.cmodel}
\title{ Generic function to summary a cluster model}
\description{ For objects of type (\code{cmodel}), this generic
summary function describes the main elements of a cluster model. }
\usage{ \method{summary}{cmodel}(object, type = "term", label = "tab:stats",
caption = NULL, ...) }
\arguments{
   \item{object}{ An object of type \code{cmodel},}
   \item{type}{ the type of xtable output in html, latex, or term,}
   \item{label}{ the label for cross reference in latex documents,}
   \item{caption}{ the caption for latex, and by extension html documents,}
   \item{...}{...}
}
\author{ Fabrice Colas }
