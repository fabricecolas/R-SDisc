\name{paste2}
\title{Concatenate strings without any space}
\description{\code{paste2} is a wrapper to paste that recursively concatenate all the string without any space nor
separators. }
\usage{paste2(...)}
\arguments{\item{...}{ one or more R objects, to be converted to character vectors.}}
\author{Fabrice Colas}
\seealso{\code{\link{paste}}}
\examples{paste2(seq(1:9))}
\alias{paste2}
