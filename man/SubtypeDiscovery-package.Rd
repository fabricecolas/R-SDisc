\name{SubtypeDiscovery-package}
\alias{SubtypeDiscovery-package}
\alias{SubtypeDiscovery}
\docType{package}
\title{Methods for Characterization, Comparison and Evaluation of Cluster
Models for Medical Research }
\description{
To compare cluster results while varying the data transformation settings or
the cluster algorithm parameters, we compute graphic- (PS) and tabular- (CSV)
reports of the essential cluster result's characteristics. We also provide
methods to cross-compare results and to evaluate the clustering stability. 

We also define complex data structures in order to compare the result of
different (1) variable transformations, to abstract on the cluster analysis (2)
settings, and to abstract on (3) clustering algorithms.

We also include two methods to evaluate the clinical relevance of clusters in a
sibling pair study (familial aggregation).
}
\details{
\tabular{ll}{
Package: \tab SubtypeDiscovery\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2008-05-08\cr
License: \tab To be defined. All rigths reserved.\cr
}



}
\author{
Fabrice Colas

Responsable : Fabrice Colas <fcolas@liacs.nl>
}
\references{
Hastie, T., Tibshirani, R., Friedman, J.: The Elements of Statistical Learning,
Data Mining, Inference, and Prediction. Springer Series in Statistics. Springer
(2001) 

Sneath, P.H.A., Sokal, R.R.: Numerical Taxonomy, The Principles and Practice of
Numerical Classification. Books in Biology. W. H. Freeman and Company (1973) 

Beyer, K., Goldstein, J., Ramakrishnan, R., Shaft, U.: When is 'nearest
neighbor' meaningful? Lecture Notes in Computer Science 1540 (1999) 217?235 

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

Eisen, M.B., Spellman, P.T., Brown, P.O., Botstein, D.: Cluster analysis and
display of genome-wide expression patterns. In: Proceedings of National Academy
of Science USA.  Volume 95. (1998) 11863-14868 

Inselberg, A.: The plane with parallel coordinates. The Visual Computer 1(2)
(1985) 69?91 Tukey, J.W.: Exploratory Data Analysis. MA: Addison-Wesley. (1977) 

Tufte, E.R.: The Visual Display of Quantitative Information. Graphics Press,
Cheshire, Connectictut Connectictut (1983) 

Tufte, E.R.: Envisioning Information. Graphics Press, Cheshire, Connectictut
(1990) 

Brewer, C.A.: 7. In: Color Use Guidelines for Mapping and Visualization.
Elsevier Science, Tarrytown, NY (1994) 123-147 

Ewens, W.J., Grant, G.R.: Statistical Methods in Bioinformatics, An
Introduction. Springer (2005) 

Gentleman, R., Carey, V.J., Huber, W., Irizarry, R.A., Dudoit, S.:
Bioinformatics and Computational Biology Solutions Using R and Bioconductor.
Springer Series in Statistics for Biology and Health. Springer (2005)
}

\keyword{ aplot }
\keyword{ dplot }
\keyword{ hplot }
\keyword{ manip }
\keyword{ classes }
\keyword{ models }
\keyword{ multivariate }
\keyword{ cluster }
\keyword{ package }

\seealso{
\code{\link[mclust:mclust-package]{mclust}} 
\code{\link[RColorBrewer:RColorBrewer-package]{RColorBrewer}} 
\code{\link[e1071:e1071-package]{e1071}} 
}
\examples{
#
# LOAD PACKAGE LIBRARY
#
library(SubtypeDiscovery)
data(setup_PD_data, setup_PD_visu_params)
# INIT 
data_4y <- analysis_PD_result <- canalysis_data <- canalysis_config <- list()
#
# LOAD SPSS DATA WHERE '1' STANDS FOR THE COLUMN WITH THE ROW NAMES (IDAA)
#
data_4y[["y1"]] <- read.delim('./data/2008-05-06_PD_data_scopayear1_n397.data', row.name=1)
data_4y[["y2"]] <- read.delim('./data/2008-05-06_PD_data_scopayear2_n346.data', row.name=1)
data_4y[["y3"]] <- read.delim('./data/2008-05-06_PD_data_scopayear3_n332.data', row.name=1)
data_4y[["y4"]] <- read.delim('./data/2008-05-06_PD_data_scopayear4_n280.data', row.name=1)
for(y in names(data_4y)){
# 
# NORMALIZE THE NAMING FOR ALL THE VARIABLES 
#
colnames(data_4y[[y]])       <- tolower(colnames(data_4y[[y]]))
data_4y[[y]]                 <- data_4y[[y]][,sort(colnames(data_4y[[y]]))]
}  
#
# SETUP PD DATA IS A MATRIX THAT DEFINES THE INDIVIDUAL VARIABLE
# TRANSFORMATIONS, PLOT COORDINATES
#
setup_PD_data[c("spesfreezing","spesaxial2"),"in_canalysis"] <- "FALSE"
for(y in names(data_4y)){
   #
   # SET THE DATA FOR EACH YEAR: INIT AND TRANSFORM DATA
   #
   #  data      := ORIGINAL DATA __MATRIX__ WITH ALL MEASUREMENTS, EVEN THOSE
   #            NOT INVOLVED IN THE CANALYSIS
   #  settings  := THE CSV-EXCEL FILE DEFINING THE VAR TRANSFORM, THE VAR FOR
   #            THE CANALYSIS, ETC. 
   #
   canalysis_data[[y]]  <- set_canalysis_data(
        data            = data_4y[[y]]          
      , settings        = setup_PD_data   
      )
   # 
   # SET THE GENERAL CONFIG FOR EACH YEAR
   # 
   #  prefix           := A PREFIX FOR ALL THE FILES GENERATED DURING THE
   #                   ANALYSIS  
   #  visu_params      := ADDITIONAL GRAPHICAL PARAMETERS
   #  stats_fun        := FUNCTIONS TO CHARACTERIZE AND EVALUTE THE CLUSTER
   #                   RESULTS
   #  bbox_threshold   := RELATIVE (%) LIMIT TO INVOLVE MODELS IN FURTHER
   #                   X-COMPARISONS
   #
   canalysis_config[[y]] <- set_canalysis_config(
      prefix=paste2("Parkinson_",y)        
      , visu_params     = setup_PD_visu_params  
      , stats_fun       = list(            
         oddratios      =          get_stats_fun(fun_name="oddratios",fun_midthreshold=mean)
         #,gen_naive_bayes =       get_stats_fun(fun_name="gen_naive_bayes",K=20)
         #,gen_knn      =          get_stats_fun(fun_name="gen_knn",K=20)
         #,gen_svm      =          get_stats_fun(fun_name="gen_svm",K=20)
         , auuc         =          get_stats_fun(fun_name="auuc"))
      , bbox_threshold  = 5                
   )
   #
   # SET THE PARAMETERS OF THE CLUSTERING ALGORITHM
   #
   calgo <- set_calgo(
      calgo = "model_based_clustering"
      , args = list(
         G_set = 3:6
         , clust_fun = mclust2 
         , model_names=c("EII", "VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV") )
   )
   #
   # DO THE ANALYSIS
   #
   analysis_PD_result[[y]] <- analysis( 
        data    = canalysis_data[[y]]
      , config  = canalysis_config[[y]]
      , calgo   = calgo) 

   # 
   # SAVE ALL MODELS AND NOT ONLY THOSE 5% LESS THAN THE BEST BIC SCORE
   #
   for(m in analysis_PD_result[[y]][["out"]])
      save_model(m, prefix=attr(canalysis_config[[y]],"prefix"))

}
rm(y)
#
# SAVE THE ?R IMAGE?
#
save.image()
}
