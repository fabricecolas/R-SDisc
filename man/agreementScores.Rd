\name{agreementScores}
\title{Agreement scores between two classifiers}
\description{Given two classification results, this function returns different statistics to assess the level of agreement: the random index, kappa, V score, the chi2 p-value and the grand total.}
\usage{agreementScores(x, y=NULL)}
\arguments{\item{x}{ classification of the first classifier}
\item{y}{ classification of the second classifier}
}
\author{Fabrice Colas}
\seealso{\code{\link{predict.SDisc}}, \code{\link{print.SDisc}}, \code{\link{summary.SDisc}}, }\keyword{distribution}
\alias{agreementScores}