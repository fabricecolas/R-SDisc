\name{evaluate_generalization}
\alias{evaluate_generalization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evaluate_generalization(data, class, fun_classifier = list(model = model_naive_bayes, predict = predict_naive_bayes), K = 7, fun_eval = stratified_traintest_split, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{class}{ ~~Describe \code{class} here~~ }
  \item{fun_classifier}{ ~~Describe \code{fun_classifier} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{fun_eval}{ ~~Describe \code{fun_eval} here~~ }
  \item{title}{ ~~Describe \code{title} here~~ }
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
function(data,class
   , fun_classifier = list(model=model_naive_bayes,predict=predict_naive_bayes)
   , K=7
   , fun_eval=stratified_traintest_split
   #, fun_eval=stratified_cv
   , title = NULL
)
{
   #
   class                   <- map(class)
   fun_transform           <- data[["fun_transform"]]
   local_data              <- cbind(data[["cc"]][,unlist(data[["sumscore_groups"]])],class=class)
   # INITIALIZATION
   indexes <- testfold     <- list()
   models <- predictions   <- contingency.tables <- perfmeasures <- list()
   G                       <- length(table(local_data$class))
   if(!(G == 0)){
      # GENERATE THE DIFFERENT STRATA FROM THE K-FOLD CROSS VALIDATION
      out_cv <- fun_eval(local_data,K)
      indexes <- out_cv[["indexes"]]
      K <- out_cv[["K"]]
      # PROCEED TO EVALUATION FOR EACH K-FOLD
      for(k in 1:K){
         # INIT THE K-TH TEST FOLD
         testfold[[k]] <- list()
         for(g in 1:G)
                  testfold[[k]] <- append(testfold[[k]],row.names(local_data[local_data$class == g,][indexes[[g]][k,],]))
         testfold[[k]] <- unlist(testfold[[k]])
         # DEFINE TRAIN AND TEST SETS
         trainset <- local_data[-pmatch(testfold[[k]], row.names(local_data)),]
         testset  <- local_data[ pmatch(testfold[[k]], row.names(local_data)),]
         #		
         for(fun_name in names(fun_transform)){
            fun_out         <- fun_transform[[fun_name]](trainset[,-match("class",dimnames(trainset)[[2]])])
            trainset        <- cbind(fun_out[["data"]][,data[["analysis_variables"]]],class=trainset$class)
            fun_out         <- fun_transform[[fun_name]](data=testset[,-match("class",dimnames(testset)[[2]])],model=fun_out[["model"]])
            testset         <- cbind(fun_out[["data"]][,data[["analysis_variables"]]],class=testset$class)
            }
         # TRAIN DIFFERENT MODELS, NAMELY NAIVE BAYES, SVM AND
         # 1-NN RK: SVM USES G(G-1) CLASSIFIERS, NAIVE BAYES AND
         # KNN ARE SINGLE MULTI-CLASS
         formula                 <- as.formula("class ~ .")
         # models[[k]]		<- fun_classifier[["model"]](formula,testset)
         models[[k]]		<- fun_classifier[["model"]](formula,trainset,testset)
         predictions[[k]]	<- fun_classifier[["predict"]](models[[k]],testset)
         contingency.tables[[k]] <- table(predictions[[k]],testset[,"class"])
         perfmeasures[[k]]	<- sum(diag(contingency.tables[[k]]))/sum(contingency.tables[[k]])
      }
      # SUMMARY STATISTICS
      perfs		<- unlist(perfmeasures)
      #
      out <- t(matrix(c(mean(perfs),1.96*sd(perfs,na.rm=TRUE),K),3,G,dimnames=list(list(title,"CI 95\%","K-folds"),as.list(1:G)),byrow=FALSE))
      out[2:G,] <- NA
      # RETURNED DATA
      return(list(models=models, predictions=predictions, contingency.tables=contingency.tables,
            perfmeasures=perfmeasures, avg=mean(perfs,na.rm=TRUE), sd=1.96*sd(perfs,na.rm=TRUE),out=out))
   }
   else
      return(NULL)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
