\name{print.SDCModel}
\title{Membership likelihood of each record into the different mixture components, along with the component mapping}
\description{Returns a matrix with \code{G}+1 columns. The \code{G} first columns represent the likelihood of each record to belong to the different mixture components. The last column corresponds the mapping of each record into one of the \code{G} components. }
\usage{print.SDCModel(x, data, ...)}
\arguments{\item{x}{ an \code{SDCModel}}
\item{data}{ an \code{SDData} container}
\item{\dots}{ (currently not used) additional parameters passed to the sub-functions}
}
\author{Fabrice Colas}
\seealso{\code{\link{modelBasedEM}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
xBestModel <- x[[bestModel(x,1)]]
print(xBestModel, data=SDData(x))}
\keyword{methods}
\keyword{print}
\alias{oddRatios}\alias{print.SDCModel}