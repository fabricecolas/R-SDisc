`evaluate_generalization` <-
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
      out <- t(matrix(c(mean(perfs),1.96*sd(perfs,na.rm=TRUE),K),3,G,dimnames=list(list(title,"CI 95%","K-folds"),as.list(1:G)),byrow=FALSE))
      out[2:G,] <- NA
      # RETURNED DATA
      return(list(models=models, predictions=predictions, contingency.tables=contingency.tables,
            perfmeasures=perfmeasures, avg=mean(perfs,na.rm=TRUE), sd=1.96*sd(perfs,na.rm=TRUE),out=out))
   }
   else
      return(NULL)
}

