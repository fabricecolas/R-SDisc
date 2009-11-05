\name{summary.bicTable}
\title{Aggregate the set of BIC scores for each combination of number of cluster, mixture model}
\description{Return an aggregated summary of the repeated data modeling performed for different combinations of number of mixture components, number of mixture model parameters and random initialization start. For each combination, it reports the best BIC score relatively to the top-ranking one along with the 95\% range of BIC scores. }
\usage{summary.bicTable(object, fun="min", bic="relativeBic", latex=FALSE, lab="bic", fmt="\%.2f", ...)}
\arguments{\item{object}{ a \code{bicTable} object}
\item{fun}{ the name of the function call for the main estimate of the table, by default the \code{min} relatie BIC is searched for}
\item{bic}{ either \code{relativeBic} -by default-, or \code{BIC}}
\item{latex}{ whether the table output should be turned into a LaTeX source code}
\item{lab}{ the label of that table in LaTeX}
\item{fmt}{ the format for the numbers in the table, by default two digits after the comma (see \code{sprintf}}
\item{\dots}{ (currently unused)}
}
\author{Fabrice Colas}
\seealso{\code{\link{bicTable}}, \code{\link{summary.bicTable}}, \code{\link{summary.SDisc}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
summary(bicTable(x))}
\keyword{cluster,methods}
\alias{summary.bicTable}