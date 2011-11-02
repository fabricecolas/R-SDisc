\name{SDCModel}
\alias{oddRatios}
\alias{oddRatiosB}
\alias{jointDistrib}
\alias{print.SDCModel}
\alias{summary.SDCModel}
\alias{SDCModel}

\title{\code{SDCModel}, \code{SDisc}'s Cluster Model class}
\description{\code{SDCModel} is the name given to \code{SDisc} Cluster Models. Subtype models from an
\code{SDisc} analysis can be accessed as elements of a list, via their name or an index integer. With \code{print},
we return the mixture component (subtype) membership likelihood per sample as well as the mapping for the most
likely mixture. The \code{summary} method characterizes the model in terms of odd ratio and joint distribution
statistics with respect to, for instance, an confounding or verificatory factor -like gender or strata- that was
not included in the subtype modeling.}

\details{
\code{print} returns a matrix with \code{G}+1 columns. The \code{G} first columns represent the likelihood of each
record to belong to the different mixture components. The last column corresponds the mapping of each record into
one of the \code{G} components. 

\code{summary} returns different statistical summaries on a  \code{SDCModel}. Among the possible summaries there
are odd ratio, chi2 statistics and joint  distribution with a categorical target variable.  Odd ratios statistics
are computed on factors of variable defined in the settings configuration file (see \code{oddGroup} in
\code{SDDataSettings}). Two types of odd ratios are implemented: the \code{oddRatios} cross-product is based on
counts of the number of occurence above and below the mean of the data distribution, while for \code{oddRatiosB}
the effect size is used in the cross-product. Concerning \code{chi2test} and \code{jointDistrib}, they involve a
comparison with a target categorical variable. 
}
\usage{
\method{print}{SDCModel}(x, data, ...)
\method{summary}{SDCModel}(object, data, type="oddRatiosB", latex=FALSE, lab="", shortStr=FALSE, ...)
}
\arguments{
\item{x}{ an \code{SDCModel}}
\item{data}{ an \code{SDData} container}
\item{object}{ an \code{SDCModel} obtained from  (see \code{SDisc} and \code{bestModel})}
\item{type}{ the type of summary to report, in: \code{"oddRatios"}, \code{"oddRatiosB"} by default, \code{"chi2test"}, \code{"jointDistrib"}}
\item{lab}{ the label of the LaTeX table}
\item{shortStr}{ whether to shorten the names of the different factors by their first letters}
\item{latex}{ whether a LaTeX formated output should be generated}
\item{\dots}{ (currently not used) additional parameters passed to the sub-functions}
}

\examples{
settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
xBestModel <- x[[bestModel(x,1)]]
print(xBestModel, data=SDData(x))
summary(x[[bestModel(x,1)]], data=SDData(x))
}
\keyword{methods}
\keyword{print}
\keyword{cluster}
\keyword{htest}

\author{Fabrice Colas}
\seealso{\code{\link{modelBasedEM}}, \code{\link{SDDataSettings}}, \code{\link{texTable}}}
