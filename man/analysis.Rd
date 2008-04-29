\name{analysis}
\alias{analysis}
\title{Generate visual-aids and table-charts in a cluster analysis with different settings}
\description{
This function acts as a wrapper to carry a cluster analysis on the provided data. 
}
\usage{
analysis(data, canalysis_variables = list(), sumscore_groups = list(),
stats_fun = list(), G_set = 1:9, modelNames = c("EII", "VII", "EEI", "EVI",
"VEI", "VVI", "EEE", "EEV", "VEV", "VVV"), formula_matrix = NULL, fun_transform
= list(identity = function(x) { return(data = I(x), model = NULL) }),
filePrefix = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Input \code{data} frame for the analysis.}
  \item{canalysis_variables}{Character vector defining the set of variables to
  consider for the cluster analysis. }
  \item{sumscore_groups}{A named list of character vectors defining the set of
  variables on which to compute odd ratios. To compute those scores, we first sum
  the variables by group into sum-scores. Then, we compare the sum-score's
  distribution in each cluster to the global distribution. To this aim, for
  each cluster, we draw 4-cells contingency table where in the rows we count
  the number of individuals member of the cluster and not, whereas in columns
  we count the number of individuals above and below the global sum-score
  average.}
  \item{stats_fun}{A named list whose elements are functions. As we estimate
  models, then we iterate over this set of functions in order to derive
  additional statistics on the cluster result. We provide each function with
  two parameters which are the \code{data} and the cluster membership
  (\code{class}). }
  \item{G_set}{A vector of small numbers which specifies the different number
  of clusters for which we should estimate models, see also \code{mclust}.}
  \item{modelNames}{A character vector that defines the types of models that we
  should estimate for our analysis, see also \code{mclust}.}
  \item{formula_matrix}{A square matrix with numbers in \{-1,0,1\} that defines
  how to adjust our data for known variability, e.g. we may want to discard the
  effect of time on our variables. So, first, the dimension names of the
  matrix should be the variable names of the \code{data}. Second, we prepare
  for each variable the regression formula to be passed to \code{lm()}, e.g.
  \code{as.formula("var_i ~ time")} with \code{var_i = 1 and time = -1}. When
  preparing the formula, +1's specify left hand side elements of the formula
  while -1's specify right hand side elements. The 0's specify variables that
  are not involved in the formula. }
  \item{fun_transform}{A named list of functions to transform the input
  \code{data}. Possible choices are \code{transform_ABSMAX, transform_AVG,
  transform_L2, transform_MEDIAN, transform_SIGMA, transform_adjust,
  transform_sibordering, transform_ALL, transform_L1, transform_MAX,
  transform_MIN , transform_addnoise, transform_longitudinal,
  transform_variables}
  }
  \item{filePrefix}{A prefix to all the output files that will be generated in
  the course of the analysis, e.g. "Parkinson" or "GARP_Study".}
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
}
\value{
  \item{orig }
  \item{cc }
  \item{analysis_outcomes }
  \item{sumscore_groups }
  \item{formula_matrix}
  \item{transformed }
  \item{fun_transform }
  \item{cluster_analysis }{
  So, the root element of \code{cluster_analysis} reports the features of the
  optimal model, namely its \code{modelName}, number of observations \code{n},
  number of dimensions \code{d}, number of clusters \code{G}, bic score
  \code{bic}, log likelihood \code{loglik}, the cluster means and the
  covariance matrix decompositions in \code{parameters}, the cluster membership
  probabilities in \code{z}, the most likely observation cluster membership in
  \code{classification}. 
  
  In \code{out}, we report \code{length(G_set)*length(modelNames)} elements
  that describe the different estimated models according to \code{mclust}. We
  overload \code{mclust} of some additional computation in order to further
  characterize the cluster results. So, each element of \code{out} reports that
  additional information within \code{stats}. Additionally, we also compute the
  footprint of each cluster result in \code{pattern}, e.g. the average-,
  median-, 5\%-, 95\%-patterns.
  
  And finally, in the root element, we also report a relative bic table, which
  we use as a basis to automatically produce a set of model cross-comparisons.}
}
\references{ ~put references to the literature/web site here ~ }
\author{Fabrice Colas}
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
#analysis(data_PD,
#   canalysis_variables  = canalysis_variables_PD,   # OUTCOMES ON WHICH THE TRANSFORMATION AND CLUSTERING SHOULD BE DONE
#   sumscore_groups      = sumscore_groups_PD,       # LOGODDS AND EVENTUALLY OTHER STATISTICS USE THIS GROUPING TO DO SUM SCORES
#   stats_fun            = stats_fun_PD,             # A LIST OF FUNCTION TO EVALUATE THE CLUSTERING RESULT
#   G_set                = 3,
#   modelNames           = c("EII","VII","EEI","EVI","VEI","VVI","EEE","EEV","VEV","VVV"),
#   fun_transform=list(
#        center=transform_AVG
#        , scale=transform_SIGMA),
#   filePrefix="Parkinson")
#}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
