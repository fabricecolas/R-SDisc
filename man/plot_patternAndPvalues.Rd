\name{plot_patternAndPvalues}
\alias{plot_patternAndPvalues}
\title{Low level plotting function for an association analysis}
\description{This function is called by the generic method
\code{plot.association_analysis} which plots the result of association analysis
performed between features and subtypes or classes. It takes essentially the
same argument than the generic method. For objects of type
(\code{association_analysis}), this generic plotting function characterizes
each subtype by its median pattern, the p-values of association between the
features and the subtype, an image summarizing the most discriminating
features. }
\usage{ plot_patternAndPvalues(a, title="", ylim=c(-2,38), ylim2=c(15,0), xlab="") }
\arguments{
   \item{a}{ An object of type \code{association_analysis},}
   \item{ylim}{ Sets the limits of the y-axis for the scores to enable comparison between the figures,}
   \item{ylim2}{ Sets the limits of the y-axis for the p-values. Notice that ylim2 is set in reverse order to get it printed at the top of each figure,}
   \item{xlab}{ Label for the x-axis,}
   \item{title}{ Text that will prefix the title of each plot.}
}
\author{ Fabrice Colas }
\seealso{\code{\link{association_analysis}}, \code{\link{plot.association_analysis}}}
