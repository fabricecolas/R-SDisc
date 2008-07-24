\name{plot.cresult}
\alias{plot.cresult}
\title{ Generic function characterizing graphically a set of cluster results}
\description{ For objects of type (\code{cresult}), this generic plotting
function draws the visual characteristicsset of a set of cluster results. This
function accepts in particular two parameters, namely \code{postscript} and
\code{query}. The first tells whether a postscript output should be generated,
whereas the second one confines the graphical characterization to a particular
queried model.}
\usage{ plot.cresult(x, device="PS", query=NULL, img_size=c(1024, 1024),...) }
\arguments{
   \item{x}{ An object of type \code{cresult},}
   \item{device}{ A character string telling the type of device to output the
   plots (in \code{"PS", "PNG"}), }
   \item{query}{ A valid reference (\code{query}) to the particular cluster
      result to illustrate. A reference can be of two types: an integer in
      \code{1:length(cresult)}, or a character string taken in
      \code{names(cresult)}, }
   \item{img_size}{ When "PNG" images are generated, default size of the plots,}
   \item{...}{ nothing }
}
\examples{
## Not run:
# LOAD DATA SET
data(wada2008)
data(wada2008_settings) 
# PREPARE CDATA
cdata_set <- set_cdata(data=wada2008, prefix="WADA2008_Sample_Analysis",
   settings=wada2008_settings)
# PREPARE THE SET OF RESULTS FOR CLUSTER MODELLING
cresult_set <- set_cresult(cdata=cdata_set,
   fun_plot=list(plot_parcoord=get_plot_fun(type="plot_parcoord"),
   plot_legend=get_plot_fun(type="plot_legend"),
   plot_image=get_plot_fun(type="plot_image"),
   plot_dendro_cluster=get_plot_fun(type="plot_dendro_cluster")),
   fun_stats=list(oddratios=get_fun_stats(fun_name="oddratios",
   fun_midthreshold=mean), auuc=get_fun_stats(fun_name="auuc")),
   nbr_top_models=5, cfun=fun_mbc_em, cfun_settings=list(modelName=c("EII",
   "VII","EEI","VEI","EVI","VVI"), G=3:6, rseed=6013:6014))
# COMPUTE THE MODEL AND ADDITIONAL STATISTICS
cresult_set <- analysis(cresult_set)
# PLOT -> PS
plot(cresult_set)
## End(Not run)
}
\author{ Fabrice Colas }
\seealso{\code{\link{set_cresult}}, \code{\link{get_plot_fun}}}
\keyword{}
