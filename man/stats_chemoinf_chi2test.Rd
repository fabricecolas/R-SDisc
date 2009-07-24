\name{stats_chemoinf_chi2test}
\alias{stats_chemoinf_chi2test}
\title{Report chi2 test statistics between a nominal attribute and the subtypes }
\description{This method is used to describe in the HTML report of a
SubtypeDiscovery analysis, the chi2 statistics between a nominal attribute and
a subtype. Values reported are the cell chi2 deviations, the sum of these
deviation and the p-value of the test.}
\usage{ stats_chemoinf_chi2test(data, class, class_col='Class', probs=seq(0,1,0.4)[2:3]) }
\arguments{
   \item{data}{ A \code{cdata} container,}
   \item{class}{ Subtype class,}
   \item{class_col}{ Nominal feature,}
   \item{probs}{ The two quantiles to help emphasize the important values in the HTML report.}
}
\author{ Fabrice Colas }
\seealso{\code{\link{stats_chemoinf_jointdistrib}}}
