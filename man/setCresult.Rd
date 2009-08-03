\name{setCresult}
\alias{setCresult}
\title{ Prepare a cluster result data structure}
\description{ This function prepares a data structure that we use to store the
set of results of a cluster analysis on a given data set. 

So, in that structure, we register one data set  (\code{cdata}) as well
as empty cluster models (\code{cmodel}) that refer to the different sketches
for the modeling.  The main element of each \code{cmodel} is actually a
prepared function call that only needs a data set as argument.}
\usage{ setCresult(cdata=NULL, cfun=clusterModelBasedEM, cfun_settings=list(modelName=c("EII",
"VII"), G=3:5, rseed=6013:6015), fun_pattern=list(mean=patternMean,
median=patternMedian, lowquant=patternLowquant, upquant=patternUpquant) ,
fun_plot=list(plotParcoord=getPlotFun(type="plotParcoord"),
plot_legend=getPlotFun(type="plot_legend"),
plot_image=getPlotFun(type="plot_image"),
plot_dendro_cluster=getPlotFun(type = "plot_dendro_cluster"),
plot_dendro_var=getPlotFun(type="plot_dendro_var")),
fun_stats=list(oddratios=getStatsFun(fun_name="oddratios")),
nTopModels=5)}
\arguments{
  \item{cdata}{ A data set of type \code{cdata}, as prepared
  by \code{setCdata}, }
  \item{cfun}{ The clustering function call,}
  \item{cfun_settings}{ A named list of vectors, which combinations are passed
  as argument to the clustering function. At this moment, as only
  \code{clusterModelBasedEM} can be passed to \code{cfun}, its default parameters are
  \code{modelName=c('EII','VII')}, \code{G=3:5} and \code{rseed=6013:6015} (Rk:
  notice that for each parameter, the length of the vector must be at least two
  for further processing of BIC patterns on a 3-dimensional array).}
  \item{fun_pattern}{ A list of functions which we use to compute different
      patterns of a cluster result,}
  \item{fun_plot}{ The different plotting function that we use to characterize
      graphically each cluster result,}
  \item{fun_stats}{ A named list containing function calls, see
  \code{getStatsFun}, }
  \item{nTopModels}{ An integer value telling how many of the most likely
      models we should consider when evaluating the clustering's consistency.  To
      help select for a particular model, we proceed to cross-comparisons of the
      cluster results. }
}
\value{
   Returns a \code{cresult} data structure of length the number of combinations
   of the \code{cfun_settings} arguments. Each element is itself a cluster
   result of type \code{cmodel}.

   For example, see demos.
}
\examples{
# LOAD THE DATA SET
data(wada2008)
# LOAD THE PRE-DEFINED SETTINGS OF WADA2008
# SEE genCdataSettings FOR MORE DETAILS
data(wada2008_settings)
wada2008_cdata <- setCdata(data=wada2008[1:100,],settings=wada2008_settings,prefix='Wada2008')
wada2008_cresult <- setCresult(cdata=wada2008_cdata)
# NOT RUN
# analysis(wada2008_cresult)
}
\author{ Fabrice Colas }
\seealso{  \code{\link{setCdata}},
\code{\link{genCdataSettings}},  \code{\link{analysis}}, .}
