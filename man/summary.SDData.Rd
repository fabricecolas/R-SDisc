\name{summary.SDData}
\title{summary.SDData}
\description{Returns a summary of the data treatments operated on the data set. For \code{mean}, \code{sd}, \code{scale}, etc, it returns the estimates. For \code{lm}, it returns the estimates of the coefficients along with their standard error, $p$-value, the $R^2$ and adjusted $R^2$ of the transformation, and the number of records on which the estimate was based on.  }
\usage{summary.SDData(object, pattern=NULL, latex=FALSE, ...)}
\arguments{\item{object}{ an data set container \code{SDData}}
\item{pattern}{ limit the summary to a subset of the data treatments (regular expression)}
\item{latex}{ whether a LaTeX formatted table output should be returned}
\item{\dots}{ (currently unused)}
}
\author{Fabrice Colas}
\seealso{\code{\link{naPattern}}, \code{\link{print.SDData}}, \code{\link{SDData}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDData(iris, settings=settings, prefix='iris')
summary(x)}
\keyword{data,htest,math,methods}
\alias{summary.SDData}