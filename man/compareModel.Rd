\name{compareModel}
\alias{compareModel}
\title{Cross compare two models}
\description{ This function takes two cluster models (\code{cmodel}) as
argument and compare them. The comparison reports the joint distribution of the
two labellings (confusion table), the level of agreement between them in terms
of the Cramer's V measure ($\chi^2$-based) but it also reports measures of
characterization (e.g log of the odds) and evaluation of the different clusters
(e.g. lambda sibs, test of goodness of fit).  }
\usage{ compareModel(x, y = NULL, fmt = c("\%d", "\%.2f")) }
\arguments{
  \item{x}{ A cluster result of type \code{cmodel},}
  \item{y}{ A cluster result of type \code{cmodel},}
  \item{fmt}{ As in \code{sprintf}.}
}
\value{ Returns a character matrix with the measurements.  }
\author{ Fabrice Colas }
\seealso{\code{\link{compareCresult}}}
