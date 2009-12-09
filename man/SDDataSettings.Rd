\name{SDDataSettings}
\title{Data settings of the SDisc analysis}
\description{Set or returns the data settings that are used to define the data: how they must be represented graphically (heatmap, parallel coordinates), what are the meaningful factors of variables, what are the data treatments. 

If a matrix is provided, it generates a 'working' settings file (settings.csv) for the current data matrix. It can be saved \code{asCSV} or returned to the standard output. 

In the case an \code{SDisc} or an \code{SDData} is provided, it returns the settings definition of the current data container or \code{SDisc} analysis.}
\usage{SDDataSettings(x, asCSV=FALSE, inCAnalysis=NULL, latex=FALSE)}
\arguments{\item{x}{ a matrix, an \code{SDData} or an \code{SDisc} object}
\item{asCSV}{ whether to save the settings  as a CSV file}
\item{inCAnalysis}{ a character vector to specify the set of variables to include in the mixture modeling}
\item{latex}{ whether a LaTeX output should be generated}}
\author{Fabrice Colas}
\seealso{\code{\link{plot.SDData}}, \code{\link{plot.SDisc}}, \code{\link{print.SDData}}, \code{\link{SDData}}, \code{\link{SDisc}}, \code{\link{summary.SDData}}, \code{\link{texTable}}, }\keyword{attribute}
\keyword{manip}
\keyword{methods}
\alias{SDDataSettings}