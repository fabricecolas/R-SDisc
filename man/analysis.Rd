\name{analysis}
\alias{analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sample workflow-method for subtype discovery}
\description{
We suggest this 'analysis' workflow as a sample for a cluster analysis. It
generates graphics reporting the result of the selected data transformation.
Then, it proceeds to a cluster analysis. Further, we characterize the cluster
results by table-reports. And we finally cross-compare the top most likely
results. 
}
\usage{
analysis(data, config, calgo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Of type 'canalysis_data' which can be constructed from the method 'set_canalyis_data()'.}
  \item{config}{Of type 'canalysis_config' which can be constructed from the method 'set_canalyis_config()'.}
  \item{calgo}{Of type 'calgo' which can be constructed from the method 'set_calgo()'.}
}
\details{
   To perform a cluster analysis, this workflow method needs three input
   parameters of class 'canalysis_data', 'canalysis_config' and 'calgo' where
   we define the way to transform the data set, the global cluster analysis
   settings and the cluster algorithm settings.

   The workflow proceeds first by an Exploratory Data Analysis (EDA) on the
   original and on the transformed data (variable_graphic_report()). We aim
   this EDA analysis to monitor the effect of our transformations. Then, we
   conduct a cluster analysis under different settings, e.g. the type of model
   or the number of clusters (calgo() with a function for main structure
   element).  We further analyze the cluster results by characterizing them
   through graphics and tables (cresult_graphic_characteristics() and
   cresult_numeric_characteristics()). Finally, as selecting models on the
   basis of BIC scores exhibiting small BIC differences, we further compare
   two-by-two the most likely models. Cross-comparisons report in particular
   the comparability of the models.

   This workflow analysis returns a regular MClust structure overloaded
   of additional computed results like graphics and tables. So, the MClust
   list-result defines its root with the elements of the bic-maximizing model.
   And the complete set of models computed in the cluster is stored in the
   'out' parameter. A 'bic_table_relative' parameter provides additional
   guidance on the relative difference in BIC score of the models, as compared
   to the optimal BIC.
}
\value{

  \item{modelName }{The name of a model in EII, VII, EEI, VEI, EVI, VVI, EEE, EEV, VEV, VVV,}
  \item{n }{The number of observations,}
  \item{d }{The number of variables in the model (or the dimension),}
  \item{G }{The number of clusters,}
  \item{bic }{The model BIC score,}
  \item{loglik }{The model log-likelihood,}
  \item{parameters }{The model parameters, i.e. principally the mean and variance estimates,}
  \item{z }{The cluster membership likelihood,}
  \item{classification }{The cluster membership,}
  \item{uncertainty }{The result of (1-cluster membership likelihood), which we can refer as the clustering 'uncertainty'.}
  \item{~ }{ }
  \item{out }{A list which stores the (MClust) models computed for the analysis under different settings,}
  \item{bic_table }{A table reporting the BIC score for all models,}
  \item{bic_table_relative }{Another BIC table reporting the relative difference of the BIC score as compared to the optimal one.}
}
\references{  
Fraley, C., Raftery, A.E.: MCLUST: Software for model-based cluster analysis.
Journal of Classi?cation 16 (1999) 297-306 

Fraley, C., Raftery, A.E.: Model-based clustering, discriminant analysis and
density estima- tion. Journal of the American Statistical Association 97 (2002)
611-631 

Fraley, C., Raftery, A.E.: Enhanced software for model-based clustering,
density estimation, and discriminant analysis: MCLUST. Journal of Classi?cation
20 (2003) 263-286 

Fraley, C., Raftery, A.E.: MCLUST version 3 for R: Normal mixture modeling and
model- based clustering. Technical Report 504, University of Washington,
Department of Statistics (September 2006) 

Banfield, J.D., Raftery, A.E.: Model-based Gaussian and non-Gaussian
clustering. Biometrics 49 (1993) 803-821 

Kass, R.E., Raftery, A.E.: Bayes factors. Journal of the American Statistical
Association 90(430) (1995) 773-? 
}
\author{ Fabrice Colas }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

## The function is currently defined as
function(data, config, calgo){
   # REDIRECT GRAPHICS -> PS
   prefix       <- attr(config,"prefix")
   # ------------
   postscript(file=paste2(prefix,"_Cluster_Analysis.ps"))
      # | VISUALIZE DISTRIB AND CORRELATION ON ORIGINAL AND TRANSFORMED DATA
      variable_graphic_report(data, which_data = "original", prefix = prefix)
      variable_graphic_report(data, prefix = prefix)
      # | CLUSTERING ANALYSIS 
      cresult           <- calgo(data, config)
      # | VISUAL CHARACTERISTICS OF THE CLUSTERS
      cresult           <- cresult_graphic_characteristics(cresult,canalysis_data = data,canalysis_config = config)
   graphics.off()
   # ------------
   # CLUSTER RESULT CHARACTERISTIZATION AND EVALUATION
   cresult              <- cresult_numeric_characteristics(cresult, canalysis_data, canalysis_config)
   # COMPARE BEST CLUSTERS  2-BY-2 AND SAVE OUTPUT INTO CSV FILES
   compare_canalysis(cresult, config=config)
   return(cresult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ models }
\keyword{ design }
\keyword{ manip }
\keyword{ hplot }
