\name{get_bic_table}
\alias{get_bic_table}
\title{ Given a set of cluster analysis results, returns a bic table}
\description{ This function returns a table that summaries the BIC scores
of a set of cluster results. The measurements reported in that table are 
the BIC scores or the scores normalized with respect to the optimal one (\%).
Finally, the format of this table is either long or wide.}
\usage{ get_bic_table(cresult, which_bic = "relative_BIC", type = "wide") }
\arguments{
  \item{cresult}{ A set of cluster results,}
  \item{which_bic}{ Which scores to return: relative to the optimal BIC score \code{relative_BIC} or simply the \code{BIC} scores,}
  \item{type}{ Which\code{type} of table format do we want, \code{long} or \code{wide}.}
}
\value{ Returns a character matrix.  }
\author{ Fabrice Colas }
\seealso{ \link[mclust:mclust-package]{bic}, \code{\link{reshape}} }
\keyword{}
