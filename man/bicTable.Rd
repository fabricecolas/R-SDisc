\name{bicTable}
\title{BIC table of an SDisc analysis}
\description{Returns the BIC table of an \code{SDisc} analysis. }
\usage{bicTable(x)}
\arguments{\item{x}{ an \code{SDisc} analysis}
}
\author{Fabrice Colas}
\seealso{\code{\link{summary.SDisc}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
bicTable(x)}
\keyword{utilities}
\alias{bicTable}