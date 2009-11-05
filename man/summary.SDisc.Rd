\name{summary.SDisc}
\title{Summarize an SDisc analysis: data report, BIC table, best models, main features of the subtypes}
\description{Returns a general summary over an \code{SDisc} analysis. Outlines the data set container characteristics, list the top-ranked models, the BIC table, and outlines oddratios of the top ranked models. }
\usage{summary.SDisc(object, q = NULL, ...)}
\arguments{\item{object}{ an \code{SDisc} analysis result}
\item{q}{ the model to summarize, e.g. one returned by \code{bestModel}}
\item{\dots}{ (currently unused) additional parameters}
}
\author{Fabrice Colas}
\seealso{\code{\link{bestModel}}, \code{\link{bicTable}}, \code{\link{print.bicTable}}, \code{\link{print.SDData}}, \code{\link{print.SDisc}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, \code{\link{summary.bicTable}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
summary(x)}
\keyword{cluster,methods}
\alias{summary.SDisc}