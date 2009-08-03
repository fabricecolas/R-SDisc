\name{stabilityAnalysis}
\alias{stabilityAnalysis}
\title{ Assess the stability to noise addition of a subtyping model}
\description{ This method assess the stability of a cluster model with respect
to noise addition. To picture the stability, we introduce in the data different
levels of noise and at each level, modeling is repeated a number of time. Next,
as models on these new data were calculated, we estimate the level of agreement
between the models. For this purpose we use the Cramer's measure of nominal
association chi2-based.

The noise levels are defined in terms of a series of numbers following:
$(1/2)^i$, with i in 1,...\code{nnoise}. The number of repeats is given by
\code{nrep}. For each repeat, the random number generator is set to a baseline
value \code{rseed} to which the index of the experiment is added, i.e. till
\code{nrep}. }
\usage{ stabilityAnalysis(x, q, rseed=6013, nnoise=10, nrep=10, ps=TRUE,
ncolors=9, minmax=c(0.5,1)) }
\arguments{
   \item{x}{ A cresult,}
   \item{q}{ Model to query for,}
   \item{rseed}{ The value to set the random number generator,}
   \item{nnoise}{ The number of noise levels to take into account,}
   \item{nrep}{ The number of replicates for each noise level,}
   \item{ps}{ Whether a postscript file must be produced,}
   \item{ncolors}{ The number of grey-colors,}
   \item{minmax}{ The minimum and maximum plot/scale values to produce comparable plots.}
}
\references{ Colas, F. et al.: Stability of Clusters for Different Time
Adjustments in Complex Disease Research. In IEEE EMBS Conference (EMBC'08),
Vancouver, Canada.  }

\author{ Fabrice Colas }

