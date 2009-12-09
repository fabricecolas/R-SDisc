\name{texTable}
\title{texTable}
\description{texTable}
\usage{texTable(x, cap="", lab="", digits=NULL, sanitize=TRUE, align = NULL, oddColor = "blueLines", longtab = FALSE, ...)}
\arguments{\item{x}{ a matrix}
\item{cap}{ the caption of the LaTeX table}
\item{lab}{ the label of the LaTeX table}
\item{digits}{ the number of digits to display}
\item{sanitize}{ whether to sanitize the row and column names values for LaTeX table (\code{xtable}). If so, define an alternative treatment function, }
\item{align}{ the alignment vector of the table}
\item{oddColor}{ the selected color to highlight alternatively the rows of the tables}
\item{longtab}{ whether a LaTeX long table environment should be used for this table}
\item{\dots}{ additional parameters to be passed to \code{xtable}}
}
\author{Fabrice Colas}
\seealso{\code{\link{print.bicTable}}, \code{\link{print.SDData}}, \code{\link{print.SDisc}}, \code{\link{SDDataSettings}}, \code{\link{summary.bicTable}}, \code{\link{summary.SDCModel}}, \code{\link{summary.SDData}}, \code{\link{summary.SDisc}}, }\keyword{manip}
\keyword{print}
\alias{texTable}