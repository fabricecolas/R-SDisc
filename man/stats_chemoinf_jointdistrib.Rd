\name{stats_chemoinf_jointdistrib}
\alias{stats_chemoinf_jointdistrib}
\title{Report joint distribution between a nominal attribute and the subtypes }
\description{This method is used to describe in the HTML report of a
SubtypeDiscovery analysis, the joint distribution between a nominal attribute
and a subtype. Values reported are the cell counts, the sum of the chi2
deviations and the p-value of the chi2 test.}
\usage{ stats_chemoinf_jointdistrib(data, class, class_col='Class') }
\arguments{
   \item{data}{ A \code{cdata} container,}
   \item{class}{ Subtype class,}
   \item{class_col}{ Nominal feature,}
}
\author{ Fabrice Colas }
\seealso{\code{\link{stats_chemoinf_chi2test}}}

