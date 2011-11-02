\name{SDStability}
\alias{SDStability}
\alias{SDStability.default}
\alias{plot.SDStability}

\title{SDStability}
\description{Evaluates the stability of a subtype model by introducing step by step more noise in the data.}
\details{

\code{plot} illustrates by color images the result of stability analyzes (\code{SDStability}).
}
\usage{
\method{SDStability}{default}(x, q, rseed=6013, nnoise=10, nrep=10)
\method{plot}{SDStability}(x, ncolors=9, minmax=c(0.5, 1), ...)
}
\arguments{\item{x}{ an \code{SDisc} analysis result}
\item{q}{ the model to assess}
\item{rseed}{ an integer to initialize the random number generator}
\item{nnoise}{ an integer to set the number of noise levels }
\item{nrep}{ an integer to indicate the number of replicates}
\item{ncolors}{ the number of colors in the grey gradient}
\item{minmax}{ the range of values for the color, by default \code{0.5} to \code{1} of the $V$ values.}
\item{\dots}{ (currently not use) additional parameters passed to the graphic functions}
}
\references{Colas, F, Meulenbelt, I, Houwing-Duistermaat, JJ, Kloppenburg, M, Watt, I, van Rooden, SM, Visser, M,
Marinus, H, van Hilten, JJ, Slagboom, PE,  Kok, JN: Stability of Clusters for Different Time Adjustments in Complex
Disease. Research 30th Annual International IEEE EMBS Conference (EMBC'08), Vancouver, Canada, 2008.
[\url{http://dx.doi.org/10.1109/IEMBS.2008.4650238}]}

\keyword{robust}
\keyword{hplot}
\keyword{methods}

\author{Fabrice Colas}
