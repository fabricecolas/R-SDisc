\name{print.bicTable}
\title{BIC table of an SDisc analysis}
\description{Returns for each combination of model, number of mixtures and random initialization seed the BIC score (Bayesian Information Criterion). Besides the absolute BIC score, we report, too, those same BIC score but relatively to absolute BIC score of the most likely model.  }
\usage{print.bicTable(x, n=NULL, modelName=NULL, G=NULL, latex=FALSE, lab="bic5", ...)}
\arguments{\item{x}{ an \code{SDisc} analysis result}
\item{n}{ either null for the whole BIC table, or a number between 0 and 1 for those models having a relative difference at most \code{n}\% worst than the most likely one (BIC), or a number greater than one returning the ordered table of the best \code{n} models. }
\item{modelName}{ limit the table to those models of type \code{modelName}}
\item{G}{ limit the table to those models with \code{G} components}
\item{latex}{ whether the table must be returned as a LaTeX source code (for dynamic reporting)}
\item{lab}{ the label of the latex table}
\item{\dots}{ additional parameters passed to \code{texTable}}
}
\author{Fabrice Colas}
\seealso{\code{\link{bicTable}}, \code{\link{SDisc}}, \code{\link{summary.SDisc}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
print(bicTable(x))}
\keyword{methods}
\keyword{print}
\alias{print.bicTable}