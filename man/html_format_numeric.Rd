\name{html_format_numeric}
\alias{html_format_numeric}
\title{Format a number into an HTML output}
\description{ To quicken the reading of the cross comparisons, we apply a
particular HTML format on the numbers. First, we color in blue and red the
numbers whose values are respectively low or high, i.e. below -1 or above +1
because these are high logodds values. Second, we remove the leading zeros in
front of the float numbers when between -1 and +1. Every cell with a zero
valued number is also replaced by an empty one. Finally, numbers are right
aligned.}
\usage{html_format_numeric(x, fmt="\%.2f", col_threshold=c(-1,1))}
\arguments{
  \item{x}{ A number: float or integer.}
  \item{fmt}{ 'sprintf'-like formating of the float numbers in the table.}
  \item{col_threshold}{ Threshold below or above which the numbers are emphasized with colors, blue below -1 and red above 1. }
}
\references{ Tufte, E.R.: The Visual Display of Quantitative Information.
Graphics Press, Cheshire, Connectictut Connectictut (1983) 

Tufte, E.R.: Envisioning Information. Graphics Press, Cheshire, Connectictut
(1990) 

Tukey, J.W.: Exploratory Data Analysis. Addison Wesley Publishing Company
(1977) }
\author{ Fabrice Colas }
\seealso{ \code{\link{print.cresult}}}
