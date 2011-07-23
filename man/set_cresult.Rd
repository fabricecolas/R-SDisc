\name{set_cresult}
\alias{set_cresult}
\title{ Prepare a cluster result data structure}
\description{ This function prepares a data structure that we use to store the
set of results of a cluster analysis on a given data set. 

So, in that structure, we register one data set  (\code{cdata}) as well
as empty cluster models (\code{cmodel}) that refer to the different sketches
for the modeling.  The main element of each \code{cmodel} is actually a
prepared function call that only needs a data set as argument.}
\usage{ set_cresult(cdata=NULL, cfun=fun_mbc_em, cfun_settings=list(modelName=c("EII",
"VII"), G=1:9, rseed=6013), fun_pattern=list(mean=patternMean,
median=patternMedian, lowquant=patternLowquant, upquant=patternUpquant) ,
fun_plot=list(plot_parcoord=get_plot_fun(type="plot_parcoord"),
plot_legend=get_plot_fun(type="plot_legend"),
plot_image=get_plot_fun(type="plot_image"),
plot_dendro_cluster=get_plot_fun(type = "plot_dendro_cluster"),
plot_dendro_var=get_plot_fun(type="plot_dendro_var")),
fun_stats=list(oddratios=get_fun_stats(fun_name="oddratios")),
fun_bic_pattern=list(mean=patternMean, median=patternMedian,
lowquant=patternLowquant, upquant=patternUpquant, sd=patternSd),
nbr_top_models=5)}
\arguments{
  \item{cdata}{ A data set of type \code{cdata}, as prepared
  by \code{set_cdata}, }
  \item{cfun}{ The clustering function call,}
  \item{cfun_settings}{ A named list of vectors, which combinations are passed
  as argument to the clustering function.}
  \item{fun_pattern}{ A list of functions which we use to compute different
      patterns of a cluster result,}
  \item{fun_plot}{ The different plotting function that we use to characterize
      graphically each cluster result,}
  \item{fun_stats}{ A named list containing function calls, see
  \code{get_fun_stats}, }
  \item{fun_bic_pattern}{ A list of functions which we use to compute different
      patterns of a cluster result,}
  \item{nbr_top_models}{ An integer value telling how many of the most likely
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
\author{ Fabrice Colas }
\seealso{  \code{\link{set_cdata}},
\code{\link{generate_cdata_settings}},  \code{\link{analysis}}, .}
\keyword{}
