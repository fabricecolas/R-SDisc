\name{summary.SDData}
\title{summary.SDData}
\description{Returns a summary of the data treatments operated on the data set. For \code{mean}, \code{sd}, \code{scale}, etc, it returns the estimates. For \code{lm}, it returns the estimates of the coefficients along with their standard error, $p$-value, the $R^2$ and adjusted $R^2$ of the transformation, and the number of records on which the estimate was based on.  }
\usage{summary.SDData(object, q=NULL, latex=FALSE,  digits = 3, ...)}
\arguments{\item{object}{ an data set container \code{SDData}}
\item{q}{ limit the summary to a subset of the data treatments (regular expression)}
\item{latex}{ whether a LaTeX formatted table output should be returned}
\item{digits}{ how many digits should be reported in the \code{SDData} summary}
\item{\dots}{ (currently unused)}
}
\author{Fabrice Colas}
\seealso{\code{\link{naPattern}}, \code{\link{print.SDData}}, \code{\link{SDData}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDData(iris, settings=settings, prefix='iris')
summary(x)}
\keyword{data}
\keyword{htest}
\keyword{math}
\keyword{methods}
\alias{summary.SDData}