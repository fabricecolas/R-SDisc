\name{model_based_clustering}
\alias{model_based_clustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
model_based_clustering(data, filePrefix = filePrefix, stats_fun, G_set = 1:9, clust_fun = mclust2, model_names = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{stats_fun}{ ~~Describe \code{stats_fun} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{clust_fun}{ ~~Describe \code{clust_fun} here~~ }
  \item{model_names}{ ~~Describe \code{model_names} here~~ }
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
function(data,filePrefix=filePrefix,stats_fun,G_set=1:9
   #, clust_fun = em_clustering  
   , clust_fun = mclust2
   , model_names=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV")
#        model_names=c("EII","VEI")
#    ,    model_names=c("EII","VVI")

        ) {
   # SUBSELECT OUTCOMES ON WHICH WE PROCEED TO THE CLUSTER ANALYSIS
   local_data                      <- data[["transformed"]][,data[["analysis_variables"]]]
   # PROCEED TO CLUSTERING
   mbc	                        <- clust_fun(data=local_data,G=G_set,modelNames=model_names)
   # FORMAT BIC TABLE
   mbc[["bic_table"]]	        <- attr(mbc[["out"]],"BIC")
   mbc[["bic_table_relative"]]     <- 100*(mbc[["bic_table"]]/max(mbc[["bic_table"]],na.rm=TRUE)-1)
   # EVENTUALLY, WRITE OUTPUT INTO CSV FILE IF FILE PREFIX IS PROVIDED
   fileNameBic                     <- paste2(filePrefix,"_Relative_BIC_Table.csv")
   if(!is.null(filePrefix))
      write.table(file=fileNameBic,mbc[["bic_table_relative"]],dec=",",sep=";")
   # SIMPLIFY OUR OWN MCLUST2 DATA STRUCTURE, ESPECIALLY "OUT" 	
   # the two fun_clust do not produce exactly the same data structure ......
   tmp_mbc_out <- attr(mbc[["out"]],"out")
   if(!is.null(tmp_mbc_out))
      mbc[["out"]] <- tmp_mbc_out
   # DERIVE STATISTICS FOR EACH CLUSTERING RESULT
   for(l in 1:length(mbc[["out"]])){
      mbc[["out"]][[l]][["stats"]] <- list()
      mbc[["out"]][[l]][["pattern"]] <- list()
      if(mbc[["out"]][[l]][["G"]] >= 2){
         if(nrow(na.omit(mbc[["out"]][[l]][["z"]])) == nrow(mbc[["out"]][[l]][["z"]])){
            class_label <- map(mbc[["out"]][[l]][["z"]])
            #
            additional_factors  <- which( !(colnames(data[["transformed"]]) \%in\% colnames(local_data)))
            mat_add_factors     <- matrix(data[["transformed"]][,additional_factors]
                                        ,nrow(data[["transformed"]])
                                        ,length(additional_factors)
                                        ,dimnames=list(
                                                row.names(data[["transformed"]])
                                                , colnames(data[["transformed"]])[additional_factors]))
            local_data_labelled <- cbind(local_data,class=class_label,mat_add_factors)
            class_count <- table(local_data_labelled$class)
            #
            mbc[["out"]][[l]][["pattern"]] <- statistical_patterns(local_data_labelled,fun_list=list(
                              avg=mean
#                                ,median=median
#                                ,s == ""lowquant=function(x){return(quantile(x,probs = seq(0, 1, 0.025))[1])}
#                                ,upquant=function(x){return(quantile(x,probs = seq(0, 1, 0.025))[40])}
                              ))
            titlePrefix <- paste2("(",mbc[["out"]][[l]][["modelName"]],",",mbc[["out"]][[l]][["G"]],") ")
            if(mbc[["out"]][[l]][["G"]] >= 3)
               graphic_characterization(mbc[["out"]][[l]][["pattern"]], 
                              class_count=class_count, 
                              titlePrefix=titlePrefix,
                              canalysis_variables = data[["analysis_variables"]])
            for(i in names(stats_fun)){
               mbc[["out"]][[l]][["stats"]][[i]] <- stats_fun[[i]](data=data,class=mbc[["out"]][[l]][["z"]])
               cat(mbc[["out"]][[l]][["G"]],mbc[["out"]][[l]][["modelName"]],", ")
            }
         }
      }
   }
   return(mbc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
