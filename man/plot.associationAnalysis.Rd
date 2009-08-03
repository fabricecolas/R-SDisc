\name{plot.associationAnalysis}
\alias{plot.associationAnalysis}
\title{ Plot association analysis between features and model}
\description{ For objects of type (\code{associationAnalysis}), this generic
plotting function characterizes each subtype by its median pattern, the
p-values of association between the features and the subtype, an image
summarizing the most discriminating features. }
\usage{ \method{plot}{associationAnalysis}(x, range_fv=NULL, p_quantile=0.95,
ylim=NA, ylim2=c(15,0), xlab="", text="", ...) }
\arguments{
   \item{x}{ An object of type \code{associationAnalysis},}
   \item{range_fv}{ The range of variables to describe graphically (case of data series),}
   \item{p_quantile}{ In the set of p-values estimated by the bootstrap method, at which level to set the threshold. For instance, setting the threshold at 0.5 would mean that 50\% of the resampled p-values are significant at 0.05 while, setting it at 0.95 means that 95\% of the resampled p-values must be significant at 0.05, which is a much stronger assumption,}
   \item{ylim}{ Sets the limits of the y-axis for the scores to enable comparison between the figures,}
   \item{ylim2}{ Sets the limits of the y-axis for the p-values. Notice that ylim2 is set in reverse order to get it printed at the top of each figure,}
   \item{xlab}{ Label for the x-axis,}
   \item{text}{ Text that will prefix the title of each plot,}
   \item{...}{ Additional arguments. }
}
\author{ Fabrice Colas }
\seealso{\code{\link{associationAnalysis}}, \code{\link{plot_patternAndPvalues}}}
