\name{modelBasedEM}
\title{Model based clustering function starting with an E-step}
\description{Proceeds to an analysis for cluster by mixture modeling (\code{mclust}). It takes as input a \code{SDCModel} (\code{x}) and a data matrix \code{data}. To init EM, a  cluster membership probability matrix (\code{z}) is drawn at random (NxG with N the number of observations and G the number of clusters). The corresponding model is estimated by an M-step and then, EM is started. }
\usage{modelBasedEM(x, data)}
\arguments{\item{x}{ an \code{SDCModel}, e.g. the first element of an \code{SDisc()[[1]]} analysis,}
\item{data}{ a data matrix, e.g. the one provided by \code{print(SDData(x)))}}
}
\references{Fraley C, Raftery AE: Model-Based Clustering, Discriminant Analysis and Density Estimation. Journal of the American Statistical Association, vol. 97, pp. 611-631, 2002. [\url{http://www.stat.washington.edu/raftery/Research/PDF/fraley2002.pdf}]Fraley C, Raftery AE, {MCLUST} Version 3 for {R}: {N}ormal Mixture Modeling and Model-Based Clustering. Technical Report 504, Department of Statistics, University of Washington, September 2006. [\url{http://www.stat.washington.edu/fraley/mclust/tr504.pdf}]}
\author{Fabrice Colas}
\seealso{\code{\link{print.SDCModel}}, \code{\link{print.SDData}}, \code{\link{SDisc}}, \code{\link{SDStability}}, }\keyword{cluster,models,multivariate}
\alias{modelBasedEM}