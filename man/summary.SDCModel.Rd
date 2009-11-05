\name{summary.SDCModel}
\title{summary.SDCModel}
\description{Returns different statistical summaries on a  \code{SDCModel}. Among the possible summaries there are odd ratio, chi2 statistics and joint  distribution with a categorical target variable. 

Odd ratios statistics are computed on factors of variable defined in the settings configuration file (see \code{oddGroup} in \code{SDDataSettings}). Two types of odd ratios are implemented: the \code{oddRatios} cross-product is based on counts of the number of occurence above and below the mean of the data distribution, while for \code{oddRatiosB} the effect size is used in the cross-product. Concerning \code{chi2test} and \code{jointDistrib}, they involve a comparison with a target categorical variable. }
\usage{summary.SDCModel(object, data, type="oddRatiosB", latex=FALSE, lab="", shortStr=FALSE, ...)}
\arguments{\item{object}{ an \code{SDCModel} obtained from  (see \code{SDisc} and \code{bestModel})}
\item{data}{ an \code{SDData} object}
\item{type}{ the type of summary to report, in: \code{"oddRatios"}, \code{"oddRatiosB"} by default, \code{"chi2test"}, \code{"jointDistrib"}}
\item{latex}{ whether a LaTeX formated output should be generated}
\item{lab}{ the label of the LaTeX table}
\item{shortStr}{ whether to shorten the names of the different factors by their first letters}
\item{\dots}{ (currently unused)}}
\author{Fabrice Colas}
\seealso{\code{\link{SDDataSettings}}, \code{\link{texTable}}, }\examples{settings <- SDDataSettings(iris)
settings['Species',] <- c(NA,FALSE, NA, NA, NA, NA)
x <- SDisc(iris, settings=settings, prefix='iris')
summary(x[[bestModel(x,1)]], data=SDData(x))}
\keyword{cluster,htest,methods}
\alias{summary.SDCModel}\alias{oddRatiosB}\alias{jointDistrib}