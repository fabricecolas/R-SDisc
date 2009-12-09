\name{print.SDisc}
\title{Assess the repeatability and/or the reproducibility of an SDisc analysis}
\description{Returns a table that cross compares the affectation of the different records into subtypes. The default settings enable to assess the repeatability of the subtypes on the two most likely models given the different combinations of number of components, mixture model, random initialization. As such, high repeatability scores assessed by the random index, V score and kappa score, indicate that there is a structure in the data that is consistently identified across the different combinations. 

In the case two different \code{SDisc} analyzes are provided in the same data, this function permits to evaluate the reproducibility of the subtypes. }
\usage{print.SDisc(x, y=NULL, m1=1, m2=2, latex=FALSE, lab="jointdistrib", ...)}
\arguments{\item{x}{ an \code{SDisc} data object}
\item{y}{ an optional second \code{SDisc} data object whose model where estimated on the same data}
\item{m1}{ rank (integer) or name (character vector) of the first model to compare. The rank is passed to \code{bestModel} to retrieve the appropriate model name. }
\item{m2}{ rank (integer) or name (character vector) of the second model to compare. The rank is passed to \code{bestModel} to retrieve the appropriate model name. }
\item{latex}{ whether the returned table must be formatted as a LaTeX table}
\item{lab}{ the label of that table}
\item{\dots}{ (currently unused)}
}
\author{Fabrice Colas}
\seealso{\code{\link{agreementScores}}, \code{\link{summary.SDisc}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
print(x)}
\keyword{methods}
\keyword{print}
\alias{print.SDisc}\alias{compareModel}