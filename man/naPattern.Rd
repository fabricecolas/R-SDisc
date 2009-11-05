\name{naPattern}
\title{Pattern of missing value of a data set}
\description{It either returns a character vector of the record ID's presenting missing values or a table (\code{latex=TRUE}) report of the missingness pattern for each record presenting at least one missing value . }
\usage{naPattern(x, latex=FALSE, ...)}
\arguments{\item{x}{ an \code{SDData} object}
\item{latex}{ by default, FALSE}
\item{\dots}{ additional arguments to be passed to \code{texTable}.}
}
\author{Fabrice Colas}
\seealso{\code{\link{SDData}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDData(iris, settings=settings, prefix='iris')
naPattern(x)}
\keyword{NA}
\alias{naPattern}