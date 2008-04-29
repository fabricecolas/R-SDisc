\name{cluster_analysis}
\alias{cluster_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cluster_analysis(data, filePrefix = NULL, G_set = 1:9, modelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV"), stats_fun = list(oddratios = function(data, class) {
    return(statistics_logodds(data, class))
}))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{modelNames}{ ~~Describe \code{modelNames} here~~ }
  \item{stats_fun}{ ~~Describe \code{stats_fun} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(data,
        filePrefix=NULL,
        G_set=1:9,
        modelNames=c("EII", "VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),
	stats_fun = list(       oddratios=function(data,class){return(statistics_logodds(data,class))}
#				,gen_naive_bayes=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_naive_bayes,predict=predict_naive_bayes),K=5))}
#				,gen_knn=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_knn,predict=predict_knn),K=5))}
#				,gen_svm=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_svm,predict=predict_svm),K=5))}
				)){
	# 
	# PROCEED TO MBC ANALYSIS ON THE TRANSFORMED DATA
	# 
	clustering <- model_based_clustering(data,filePrefix=filePrefix,stats_fun=stats_fun,G_set=G_set)	
	return(clustering)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
