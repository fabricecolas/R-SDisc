\name{print.SDData}
\title{Matrix used to perform the cluster analysis}
\description{Returns the data matrix that results from the initialization (e.g. complete cases) and  the data treatments applied on the different variables, as defined in the data settings configuration file. 

To verify the data treatments, an \code{rseed} can be provided to place in parallel a random extract of the data matrix before and after the data treatments. The \code{range} parameter gives the number of rows and columns to extract randomly. }
\usage{print.SDData(x, rseed=NULL, range=1:3, latex=FALSE, ...)}
\arguments{\item{x}{ an \code{SDData} object}
\item{rseed}{ an integer to set the random number generator that will select randomly a set of rows and columns to see before and after the data treatment, the data matrix}
\item{range}{ a sequence of integers to subset the randomly ordered column and row names vectors}
\item{latex}{ whether the table  must be returned into a LaTeX code }
\item{\dots}{ (currently unused) additional parameters to be passed to the subfunctions}
}
\author{Fabrice Colas}
\seealso{\code{\link{modelBasedEM}}, \code{\link{SDData}}, \code{\link{SDDataSettings}}, \code{\link{summary.SDData}}, }\keyword{methods,print}
\alias{print.SDData}