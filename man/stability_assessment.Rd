\name{stability_assessment}
\alias{stability_assessment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stability_assessment(data = na.omit(as.data.frame(df)), analysis_variables = NULL, sumscore_groups = NULL, stats_fun = NULL, fun_transform = list(center = transform_AVG, scale = transform_SIGMA), G_set = 1:9, query_model = "VVI,5", modelNames = c("EII", "VII", "EEI", "EVI", "VEI", "VVI", "EEE", "EEV", "VEV", "VVV"), noise_vector = (1/2^(1:10)), K = 10, filePrefix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{analysis_variables}{ ~~Describe \code{analysis_variables} here~~ }
  \item{sumscore_groups}{ ~~Describe \code{sumscore_groups} here~~ }
  \item{stats_fun}{ ~~Describe \code{stats_fun} here~~ }
  \item{fun_transform}{ ~~Describe \code{fun_transform} here~~ }
  \item{G_set}{ ~~Describe \code{G_set} here~~ }
  \item{query_model}{ ~~Describe \code{query_model} here~~ }
  \item{modelNames}{ ~~Describe \code{modelNames} here~~ }
  \item{noise_vector}{ ~~Describe \code{noise_vector} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
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
function(
     data               = na.omit(as.data.frame(df))
   , analysis_variables  = NULL
   , sumscore_groups    = NULL
   , stats_fun          = NULL
   , fun_transform      = list(
           center       = transform_AVG
         , scale        = transform_SIGMA)   
   , G_set              = 1:9
   , query_model        = "VVI,5"
   , modelNames         = c("EII","VII","EEI","EVI","VEI","VVI","EEE","EEV","VEV","VVV")
   , noise_vector       = (1/2^(1:10))
   , K                  = 10
   , filePrefix         = NULL 
   ){
   #
   canalysis_results               <- list()
   for(adj in names(adj_params)){
      canalysis_results[[adj]]     <- list()
      for(noise_idx in 1:length(noise_vector)){
         canalysis_results[[adj]][[noise_idx]] <- list()
         # ADD A NOISE FUNCTION TO THE LIST OF DATA TRANSFORMATIONS
         fun_transform_set      <- adj_params[[adj]]$transform
         mat_formula            <- adj_params[[adj]]$formula_matrix
         fun_transform[["time_adj"]] <-  function(data, model = NULL, type  = "lm"){ return(
                                                transform_adjust(data, model, type, f_matrix = mat_formula, transform=fun_transform_set))}
         fun_transform[["center"]]<- transform_AVG
         fun_transform[["scale"]] <- transform_SIGMA
         fun_transform[["noise"]] <- function(data, model= NULL, type="addnoise"){return(
                                                transform_addnoise(data, model, canalysis_variables = analysis_variables, type, 
                                                relative_noise_amount = noise_vector[noise_idx], rand_seed = 6013+3*i))}
         filePrefix_local <- paste2(filePrefix,noise_vector[noise_idx],"_",adj,"_")
         # PROCEED TO CLUSTER ANALYSIS 10 TIMES WITH A DIFFERENT SEED
         for(i in 1:10)
            canalysis_results[[adj]][[noise_idx]][[i]] <- analysis(
                                                              data                   = data
                                                            , analysis_variables      = analysis_variables 
                                                            , sumscore_groups        = sumscore_groups
                                                            , stats_fun              = stats_fun
                                                            , fun_transform          = fun_transform
                                                            , G_set                  = G_set
                                                            , K                      = K
                                                            , filePrefix             = filePrefix_local )
         # PROCEED TO CORRELATION EVALUATION
         b <- list(auuc=c(),cramerv=c())
         for(i in 1:10){
            for(j in 1:10){
               if(j > i){
                  a <- compare_canalysis(
                           canalysis_results[[adj]][[noise_idx]][[i]],
                           canalysis_results[[adj]][[noise_idx]][[j]], 
                           query=query_model, filePrefix=NULL)
                  if(length(a) > 0){
                     a <- a[[1]]
                     b[["cramerv"]] <- c(b[["cramerv"]],a["Cramer's V",1])
                  }
               }
            }
         }
         # PRESERVE THE RESULT
         canalysis_results[[adj]][[noise_idx]][["result"]] <- b
         canalysis_results[[adj]][[noise_idx]][["result"]][["mean_cramerv"]] <- mean(as.numeric(b[["cramerv"]]),na.rm=TRUE)
         canalysis_results[[adj]][[noise_idx]][["result"]][["sd_cramerv"]] <- sd(as.numeric(b[["cramerv"]]),na.rm=TRUE)
      }
   }
   #
   # QUANTILE'S DATA FRAME
   #
   # DF
   quantile_df <- matrix(NA,length(noise_vector),11,dimnames=list(100*noise_vector,c("0\%","10\%","20\%","30\%","40\%","50\%","60\%","70\%","80\%","90\%","100\%")))
   # ARRAY
   quantile_a <- quantile_df
   for(adj in 2:length(adj_params))
      quantile_a <- abind(quantile_a,quantile_df,along=3)
   dimnames(quantile_a)[[3]] <- names(adj_params)
   # FILL ARRAY VALUES 
   for(adj in names(adj_params))
      for(noise_idx in 1:length(noise_vector))
         quantile_a[noise_idx,,adj] <- quantile(as.numeric(canalysis_results[[adj]][[noise_idx]][["result"]]$cramerv),probs=seq(0,1,0.1))
   #
   ncolors <- 9
   local_breaks <- seq(min(quantile_a),max(quantile_a),(max(quantile_a)-min(quantile_a))/ncolors)
   local_colors <- brewer.pal(ncolors,"Greys")
   postscript(file=paste2(filePrefix,"_stability.ps"),horizontal=FALSE,paper="a4",pagecentre=TRUE)
   par(mfrow=c(2,2))
   for(adj in names(adj_params)){
      image(quantile_a[,,adj],col=local_colors,breaks=local_breaks,axes=FALSE,ylab=adj,cex.lab=1.7)
      contour(quantile_a[,,adj],add=TRUE,labcex=1.3)
      axis(1,at=seq(0,1,0.11),labels=sprintf("\%.1f\%\%",as.numeric(row.names(quantile_a))),las=2,cex.axis=1.7)
      axis(2,at=seq(0,1,0.1),labels=colnames(quantile_a),cex.axis=1.7)
      }
   graphics.off()
   return(list(canalysis=canalysis_results,quantile_a))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
