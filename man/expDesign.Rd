\name{expDesign}
\title{Charts the Joint Distribution of Groups of Factors Two by Two}
\description{The method enables to report the joint distributions of two or more factors with each others.  If
\code{gr1} and \code{gr2} are both \code{NULL}, then each column of \code{x} is compared two-by-two. If \code{gr1}
is not \code{NULL} but \code{gr2} is, then the different factors of \code{gr1} are compared one to another.
Finally, if \code{gr2} is not null, the factors of \code{gr1} are compared with those of \code{gr2}.  }
\usage{expDesign(x, gr1 = NULL, gr2 = NULL, tex = FALSE, exclude = NULL)}
\arguments{
\item{x}{a data frame or matrix, }
\item{gr1}{a set of factor names from \code{x}, \code{NULL} by default}
\item{gr2}{a second set of factor names from \code{x} to compare against the factors from \code{gr1}, \code{NULL}
by default}
\item{tex}{if \code{TRUE} a (Sweave-friendly) latex output is generated; if \code{FALSE} tables are directly
returned to the standard output}
\item{exclude}{should missing values be considered as level; by default yes \code{NULL} -see help of \code{table}-.
} }
\author{Fabrice Colas}
\examples{
expDesign(OrchardSprays, gr1='treatment', gr2='rowpos')
expDesign(OrchardSprays, gr1='treatment', gr2='rowpos', tex=TRUE)
}
\alias{expDesign}
