\name{bestModel}
\title{Most likely mixture models of an SDisc analysis}
\description{Taking as input an \code{SDisc} analysis, returns a  character vector of the most likely models ordered by decreasing likelihood of length \code{n}. Optional parameters \code{modelName} and \code{G} may act as filter to subset the most likely model given those parameters. }
\usage{bestModel(x, n=NULL, modelName=NULL, G=NULL)}
\arguments{\item{x}{ an \code{SDisc} object}
\item{n}{ the number of model to return, by default 5 as specified in \code{SDisc}}
\item{modelName}{ filter the returned models matching \code{modelName}}
\item{G}{ filter the returned models to those matching with \code{G}}
}
\author{Fabrice Colas}
\seealso{\code{\link{summary.SDCModel}}, }\keyword{arith,character}
\alias{bestModel}