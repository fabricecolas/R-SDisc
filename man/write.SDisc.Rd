\name{write.SDisc}
\title{write.SDisc}
\description{For the 5 best models (by default, see \code{nTopModel} in \code{SDisc}), writes as CSV files the cluster membership matrix of the different records along with their affectation into the different clusters (see \code{print.SDCModel}). By default, \code{SDisc} performs this task and save the models. }
\usage{write.SDisc(x, q=NULL)}
\arguments{\item{x}{ an \code{SDisc} analysis result}
\item{q}{ a NULL, numeric or character vector as provided by \code{bestModel}.}
}
\author{Fabrice Colas}
\seealso{\code{\link{print.SDisc}}, \code{\link{SDisc}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
write.SDisc(x)}
\keyword{cluster}
\keyword{print}
\keyword{utilities}
\alias{write.SDisc}