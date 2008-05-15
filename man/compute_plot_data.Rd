\name{compute_plot_data}
\alias{compute_plot_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Compute plot-data patterns like average, median or quantiles}
\description{ Given the function provided as parameter 'fun', it computes the
patterns for each cluster. For instance, providing mean as parameter would give
the average of each cluster, i.e. its center. 

The function returns a matrix of dimension dxG, with d the number of variables
in the data set and G the number of clusters. 

Remark that \code{compute_plot_data()} does not restrict to the variables for
the cluster analysis but to all of them and we treat the missing values by
omitting them.  }
\usage{ compute_plot_data(data, fun = mean) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A data matrix with a \code{class} column, referring to the cluster membership,}
  \item{fun}{ A function like \code{mean, sd, median} or even
      {function(x){return(quantile(x,probs=0.025,na.rm=TRUE))}} for the 2.5\%
      quantile-pattern of the cluster.}
}
\details{ }
\value{
  \item{comp1 }{Description of 'comp1'}
}
\references{ }
\author{ Fabrice Colas }
\note{ ~Make other sections like Warning with \section{Warning }{....} ~ }
\seealso{ \code{\link{set_plot_data}}, \code{\link{get_coloring_scheme}}}
\examples{ }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
