\name{predict.SDisc}
\title{Estimates a set of subtypes (mixture models) on new data}
\description{Retrieve the most likely mixture model of an \code{SDisc} analysis and use it for prediction on new data. Enables to conduct a validation analysis on the profiles obtained from a first \code{SDisc} analysis with the ones from a second \code{SDisc} analysis. To be used in conjunction with \code{print.SDisc} that can cross-compare the profile affectation of two models on the same data. }
\usage{predict.SDisc(object, newdata, ...)}
\arguments{\item{object}{ an \code{SDisc} object}
\item{newdata}{ an \code{SDData} or \code{SDisc} object}
\item{\dots}{ (currently not used) additional parameters to be passed to the subfunctions}
}
\author{Fabrice Colas}
\seealso{\code{\link{print.SDisc}}, \code{\link{SDData}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, \code{\link{summary.SDisc}}, }\keyword{classif}
\keyword{cluster}
\keyword{methods}
\keyword{models}
\keyword{multivariate}
\alias{predict.SDisc}