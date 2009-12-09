\name{plot.SDisc}
\title{Graphical report of an SDisc analysis}
\description{Draws a graphical report of an \code{SDisc} analysis. By default (\code{q=NULL}), the best model is represented by different plots \code{type} whose configuration is defined in \code{SDDataSettings(x)}. Thus, it will draw the parallel coordinate, legend, image, dendrogram on the subtypes/profiles and another on the variables. Generated graphics are made into PDF's. If \code{latex=TRUE}, the LaTeX code to enclose the graphics into a dynamic reporting tool is generated. Additional parameters  enable to influence the color scales, range of values and of features.}
\usage{plot.SDisc(x, q=NULL, type=c("plotParcoord", "plotLegend", "plotImage", "plotDendroCluster", "plotDendroVar"), latex=FALSE, title=NULL, xlim=c(-3, 3), zlim=c(-2, 2), xy=c(-2.2, 0), pattern=mean, cex=0.7, colGrad=rev(brewer.pal(9, "RdBu")), rangeFV=NULL, lab=NULL, ...)}
\arguments{\item{x}{ an \code{SDisc} analysis}
\item{q}{ an numeric value telling how many top ranking models to characterize graphically (\code{bestModel}, a character vector referring to the names of an \code{SDCModel}, by default set to NULL that is, the 5 most likely models)}
\item{type}{ a character vector in \code{'plotParcoord', 'plotLegend', 'plotImage', 'plotDendroCluster', 'plotDendroVar'}}
\item{latex}{ either TRUE or FALSE, whether the LaTeX code must be reported on the standard output for dynamic report generation (\code{Sweave})}
\item{title}{ the title of the graphics}
\item{xlim}{ the x-limits of the parallel coordinate plots}
\item{zlim}{ the z-limits for the color gradient in the image}
\item{xy}{ the xy-location of the legend}
\item{pattern}{ the name of the function to calculate the characteristic pattern, by default \code{mean}}
\item{cex}{ the character expansion numeric value, by default \code{0.7} 70\%}
\item{colGrad}{ a character vector of the colors to use in the color image}
\item{rangeFV}{ the range of features when plotting series of values}
\item{lab}{ the label of the LaTeX figure, \code{NULL} by default}
\item{\dots}{ (currently not used) additional graphical parameters passed to the plotting functions}
}
\author{Fabrice Colas}
\seealso{\code{\link{plot.SDisc}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
#plot(x)}
\keyword{hplot}
\keyword{methods}
\alias{plot.SDisc}\alias{plotSeries}\alias{plotParcoord}\alias{plotDendroCluster}\alias{plotDendroVar}\alias{plotImage}\alias{plotLegend}