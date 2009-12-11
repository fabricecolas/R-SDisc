\name{SDisc}
\title{Perform an SDisc analysis}
\description{Performs and returns an \code{SDisc} analysis on the \code{data}. This analysis involves repeated mixture modeling for different combinations of number of components, number of mixture model parameters and random initialization start. The \code{SDisc} results contains a data set container (\code{SDData}) which stores the original data and which may, e.g., limit the cluster analysis to a few variables (see \code{SDData}). It contains, too, the parameters of the different models estimated and a BIC table summarizing their likelihood and rank. Generic plot, print and summary function enable to visualize and summarize the results (see \code{plot.SDisc}, \code{print.SDisc} and \code{summary.SDisc}).  }
\usage{SDisc(data, cfun="modelBasedEM", cFunSettings=list(modelName = c("EII", "VII"), G = 3:5, rseed = 6013:6015), nTopModels=5, nnodes=1, ...)}
\arguments{\item{data}{ a data matrix (with its \code{settings} file) or an \code{SDData} instance}
\item{cfun}{ the name of the cluster algorithm}
\item{cFunSettings}{ the set of parameters of the cluster algorithm}
\item{nTopModels}{ the number of top-ranking models} 
\item{nnodes}{ the number of nodes in the case of parallel computing}
\item{\dots}{ additional parameters passed to the \code{SDData} function when a \code{data} matrix is provided. Such parameters may be \code{settings} and \code{prefix}.}}
\references{Colas, F, Meulenbelt, I, Houwing-Duistermaat, JJ, Kloppenburg, M, Watt, I, van Rooden, SM, Visser, M, Marinus, H, Cannon, EO, Bender, A, van Hilten, JJ, Slagboom, PE, Kok, JN: A Scenario Implementation in R for Subtype Discovery Examplified on Chemoinformatics Data. Leveraging Applications of Formal Methods, Verification and Validation (ISoLA'08), October 13-15, 2008. [\url{http://dx.doi.org/10.1007/978-3-540-88479-8_48}]
\author{Fabrice Colas}
\seealso{\code{\link{plot.SDisc}}, \code{\link{plot.SDStability}}, \code{\link{predict.SDisc}}, \code{\link{print.SDisc}}, \code{\link{SDData}}, \code{\link{SDDataSettings}}, \code{\link{SDisc}}, \code{\link{summary.SDisc}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
#plot(x)
print(x)
summary(x)}
\keyword{cluster}
\keyword{math}
\keyword{methods}
\keyword{multivariate}
\alias{SDisc}