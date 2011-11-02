\name{SDisc}
\alias{SDisc}
\alias{SDisc.default}
\alias{compareModel}
\alias{plot.SDisc}
\alias{plotSeries}
\alias{plotParcoord}
\alias{plotDendroCluster}
\alias{plotDendroVar}
\alias{plotImage}
\alias{plotLegend}
\alias{predict.SDisc}
\alias{print.SDisc}
\alias{summary.SDisc}
\title{SDisc to discover homogeneous subtypes in data}
\description{Performs and returns an \code{SDisc} analysis on the \code{data}. This analysis involves repeated mixture modeling for different combinations of number of components, number of mixture model parameters and random initialization start. The \code{SDisc} results contains a data set container (\code{SDData}) which stores the original data and which may, e.g., limit the cluster analysis to a few variables (see \code{SDData}). It contains, too, the parameters of the different models estimated and a BIC table summarizing their likelihood and rank. Generic plot, print and summary function enable to visualize and summarize the results (see \code{plot.SDisc}, \code{print.SDisc} and \code{summary.SDisc}).  }
\usage{
\method{SDisc}{default}(x, cfun="modelBasedEM", cFunSettings=list(modelName = c("EII", "VII"), G = 3:5, rseed = 6013:6015), nTopModels=5, nnodes=1, ...)


\method{plot}{SDisc}(x, q=NULL, type=c('plotParcoord', 'plotLegend', 'plotPC1', 'plotPC2', 'plotDendroCluster',
   'plotDendroVar'), latex=FALSE, title=NULL, xlim=c(-3, 3), zlim=c(-2, 2), xy=c(-2.2, 0), pattern=mean, cex=0.7,
   colGrad=rev(brewer.pal(9, "RdBu")), rangeFV=NULL, lab=NULL, ...)
\method{predict}{SDisc}(object, newdata, ...)
\method{print}{SDisc}(x, y=NULL, m1=1, m2=2, latex=FALSE, lab="jointdistrib", ...)
\method{summary}{SDisc}(object, q = 1, ...)
}
\arguments{
\item{x}{ a data matrix (with its \code{settings} file), an \code{SDData} instance or an \code{SDisc} data object}
\item{cfun}{ the name of the cluster algorithm}
\item{cFunSettings}{ the set of parameters of the cluster algorithm}
\item{nTopModels}{ the number of top-ranking models} 
\item{nnodes}{ the number of nodes in the case of parallel computing}
\item{object}{ an \code{SDisc} analysis result}
\item{q}{ an numeric value telling how many top ranking models to characterize graphically (\code{bestModel}, a
character vector referring to the names of an \code{SDCModel}, by default set to NULL that is, the 5 most likely
models)}
\item{y}{ an optional second \code{SDisc} data object whose model where estimated on the same data}
\item{m1}{ rank (integer) or name (character vector) of the first model to compare. The rank is passed to \code{bestModel} to retrieve the appropriate model name. }
\item{m2}{ rank (integer) or name (character vector) of the second model to compare. The rank is passed to \code{bestModel} to retrieve the appropriate model name. }
\item{latex}{ either TRUE or FALSE, whether the LaTeX code must be reported on the standard output for dynamic report generation (\code{Sweave})}
\item{lab}{ the label of that table}
\item{type}{ a character vector in \code{'plotParcoord', 'plotLegend', 'plotImage', 'plotDendroCluster', 'plotDendroVar'}}
\item{title}{ the title of the graphics}
\item{xlim}{ the x-limits of the parallel coordinate plots}
\item{zlim}{ the z-limits for the color gradient in the image}
\item{xy}{ the xy-location of the legend}
\item{pattern}{ the name of the function to calculate the characteristic pattern, by default \code{mean}}
\item{cex}{ the character expansion numeric value, by default \code{0.7} 70\%}
\item{colGrad}{ a character vector of the colors to use in the color image}
\item{rangeFV}{ the range of features when plotting series of values}
\item{newdata}{ an \code{SDData} or \code{SDisc} object}
\item{\dots}{ additional parameters passed to the \code{SDData} function when a \code{data} matrix is provided.
Such parameters may be \code{settings} and \code{prefix}.}
}
\details{

\code{plot} characterizes graphically the best model (\code{q=NULL}) of an \code{SDisc} analysis through parallel
coordinate, legend, image, sample- and variable-dendrograms. Some plots like heatmap and parallel coordinates admit
parameterization, which can be specified in the data description (see \code{SDDataSettings(x)}). Generated graphics
are in PDF and Sweave output can be switched on with \code{latex=TRUE}. Additional parameters  influence the color
scales, range of values.

\code{predict} estimates a set of subtypes (mixture models) on new data based on the most likely mixture model from
an \code{SDisc} analysis. Used in combination with \code{print},  the new data can help validate a previous model
of subtypes. 

\code{print} reports the joint distribution between pairs of subtyping results, the agreement (kappa, V
chi2-based measure, random index). It is also possible to compare affectation between \code{SDisc} analyses
with different parameters but with the same data container. In which case, it illustrates how, under different
settings, the subtyping analyses agree. 

\code{summary} lists the dataset characteristics, the top-ranked models, the BIC table, and oddratios statistics
for the top ranked models.
}
\references{Colas, F, Meulenbelt, I, Houwing-Duistermaat, JJ, Kloppenburg, M, Watt, I, van Rooden, SM, Visser, M,
Marinus, H, Cannon, EO, Bender, A, van Hilten, JJ, Slagboom, PE, Kok, JN: A Scenario Implementation in R for
Subtype Discovery Examplified on Chemoinformatics Data. Leveraging Applications of Formal Methods, Verification and
Validation (ISoLA'08), October 13-15, 2008. [\url{http://dx.doi.org/10.1007/978-3-540-88479-8_48}]

Colas, F, Meulenbelt, I, Houwing-Duistermaat, JJ, Kloppenburg, M, Watt, I, van Rooden, SM, Visser, M, Marinus, H,
van Hilten, JJ, Slagboom, PE,  Kok, JN: Stability of Clusters for Different Time Adjustments in Complex Disease.
Research 30th Annual International IEEE EMBS Conference (EMBC'08), Vancouver, Canada, 2008.
[\url{http://dx.doi.org/10.1109/IEMBS.2008.4650238}]

Fraley C, Raftery AE: Model-Based Clustering, Discriminant Analysis and Density Estimation. Journal of the American
Statistical Association, vol. 97, pp. 611-631, 2002.
[\url{http://www.stat.washington.edu/raftery/Research/PDF/fraley2002.pdf}]

Fraley C, Raftery AE, {MCLUST} Version 3 for {R}: {N}ormal Mixture Modeling and Model-Based Clustering. Technical
Report 504, Department of Statistics, University of Washington, September 2006.
[\url{http://www.stat.washington.edu/fraley/mclust/tr504.pdf}]
}
\examples{
settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
### do not run
#plot(x)
print(x)
summary(x)
}

\keyword{classif}
\keyword{cluster}
\keyword{math}
\keyword{methods}
\keyword{models}
\keyword{multivariate}
\keyword{print}
\keyword{hplot}

\seealso{\code{\link{SDStability}}, \code{\link{SDData}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, 
\code{\link{bestModel}}, \code{\link{bicTable}}, \code{\link{texTable}}, \code{\link{agreementScores}}}

\author{Fabrice Colas}
