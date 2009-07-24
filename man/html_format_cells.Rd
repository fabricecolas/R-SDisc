\name{html_format_cells}
\alias{html_format_cells}
\title{Quantile-based html number formatting}
\description{ To quicken the reading of the cross comparisons, we apply a
particular HTML format on the numbers. First, we color in blue and red the
numbers whose values are respectively low or high, i.e. below -1 or above +1
because these are high logodds values. Second, we remove the leading zeros in
front of the float numbers when between -1 and +1. Every cell with a zero
valued number is also replaced by an empty one. Finally, numbers are right
aligned.}
\usage{ html_format_cells(x, fmt="\%.2f", quant_threshold=quantile(c(-1,1),probs=seq(0,1,0.333))) }
\arguments{
  \item{x}{ A number: float or integer.}
  \item{fmt}{ 'sprintf'-like formating of the float numbers in the table.}
   \item{quant_threshold}{ ,}
   \item{probs}{ ,}
}
\references{ Tufte, E.R.: The Visual Display of Quantitative Information.
Graphics Press, Cheshire, Connectictut Connectictut (1983) 

Tufte, E.R.: Envisioning Information. Graphics Press, Cheshire, Connectictut
(1990) 

Tukey, J.W.: Exploratory Data Analysis. Addison Wesley Publishing Company
(1977) }
\author{ Fabrice Colas }
\seealso{\code{\link{html_format_cells}}}
