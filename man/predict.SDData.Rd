\name{predict.SDData}
\title{Transform a new data set based on the transformations from another SDData}
\description{Given a previously treated dataset, transform new data using the same estimates.  }
\usage{predict.SDData(object, newdata, prefix = "Newdata", subset = NULL, ...)}
\arguments{\item{object}{ the previously calculated \code{SDData} object}
\item{newdata}{ a new dataset on which to apply the same transformation estimated on the first \code{SDData}}
\item{prefix}{ the prefix of the new \code{SDData}}
\item{subset}{ a row index to subset the records of the dataset}
\item{\dots}{ (currently unused) additional parameters to be passed to the subfunctions}}
\author{Fabrice Colas}
\seealso{\code{\link{predict.SDisc}}, \code{\link{SDData}} }
\keyword{manip}
\keyword{multivariate}
\alias{predict.SDData}
