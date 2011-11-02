\name{bicTable}
\alias{summary.bicTable}
\alias{bicTable}
\alias{bicTable.default}
\alias{print.bicTable}
\title{BIC table of an SDisc analysis}
\description{Returns for each combination of model, number of mixtures and random initialization seed the BIC score
(Bayesian Information Criterion). Besides the absolute BIC score, we report, too, those same BIC score but
relatively to absolute BIC score of the most likely model. }
\details{
\code{summary} returns an aggregated summary of the repeated data modeling performed for different combinations of
number of mixture components, number of mixture model parameters and random initialization start. For each
combination, it reports the best BIC score relatively to the top-ranking one along with the 95\% range of BIC
scores.
}
\usage{
\method{bicTable}{default}(x, ...)
\method{print}{bicTable}(x, n=NULL, modelName=NULL, G=NULL, latex=FALSE, lab="bic5", ...)
\method{summary}{bicTable}(object, fun="min", bic="relativeBic", latex=FALSE, lab="bic", fmt="\%.2f", ...)
}
\arguments{
\item{x}{ an \code{SDisc} analysis}
\item{n}{ either null for the whole BIC table, or a number between 0 and 1 for those models having a relative difference at most \code{n}\% worst than the most likely one (BIC), or a number greater than one returning the ordered table of the best \code{n} models. }
\item{modelName}{ limit the table to those models of type \code{modelName}}
\item{G}{ limit the table to those models with \code{G} components}
\item{latex}{ whether the table must be returned as a LaTeX source code (for dynamic reporting)}
\item{lab}{ the label of the latex table}
\item{object}{ a \code{bicTable} object}
\item{fun}{ the name of the function call for the main estimate of the table, by default the \code{min} relatie BIC is searched for}
\item{bic}{ either \code{relativeBic} -by default-, or \code{BIC}}
\item{fmt}{ the format for the numbers in the table, by default two digits after the comma (see \code{sprintf})}
\item{\dots}{ additional parameters passed to subfunctions, e.g. \code{texTable}}
}
\examples{
settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
bicTable(x)
}

\keyword{cluster}
\keyword{methods}
\keyword{print}
\keyword{utilities}

\seealso{\code{\link{bicTable}}, \code{\link{SDisc}}, \code{\link{texTable}}}

\author{Fabrice Colas}
