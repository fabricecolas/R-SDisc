\name{plot.SDData}
\title{Graphical report of an SDData data set container}
\description{Reports as PDF the boxplots and histograms for each variable of the data set container. If \code{latex} is set to TRUE, then the LaTeX code to include the different figures is returned into the standard output (for use with \code{Sweave}). }
\usage{plot.SDData(x, q=NULL, est = 1, zlim = c(-2, 2), latex=FALSE, ...)}
\arguments{\item{x}{ an \code{SDData} object}
\item{q}{ limit the summary to a subset of the data treatments (regular expression). When a character vector is provided, the data treatments matching the regular expression are plotted side by side in an image.}
\item{est}{ TODO}
\item{zlim}{ TODO}
\item{latex}{ either TRUE or FALSE, by default FALSE and simply generates the PDF of the figures}
\item{\dots}{ (currently not used) additional parameters to be passed to the plotting functions.}
}
\author{Fabrice Colas}
\seealso{\code{\link{print.SDData}}, \code{\link{SDData}}, \code{\link{summary.SDData}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDData(iris, settings=settings, prefix='iris')
#plot(x)}
\keyword{hplot}
\keyword{methods}
\alias{plot.SDData}