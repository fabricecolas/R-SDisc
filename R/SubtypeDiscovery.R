`analysis` <-
function(cresult, doPSOutput = TRUE, csv = TRUE, html = TRUE){
   plot(attr(cresult,"cdata"), doPSOutput = doPSOutput)
   cresult      <- fun_cmodel(cresult)
   plot(cresult, doPSOutput = doPSOutput)
   print(cresult,html = TRUE)
   return(cresult)
}

`stats_auuc` <-
function(class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}

`compare_cresult` <-
function(cresult,cresult2=NULL,query=NULL){
   prefix               <- attr(cresult,"prefix")
   nbr_top_models       <- attr(cresult,"nbr_top_models")
   cross_comparison_results <- list()
   best_models          <- get_best_models(cresult)
   model_pairs          <- t(combn(best_models,m=2))
   if(nrow(model_pairs)>=1){
      cnames            <- gsub(",","_",paste("(",apply(model_pairs,1,paste,collapse="),("),")",sep=""))
      for(mp_idx in 1:nrow(model_pairs)){
         x              <- cresult[[which(model_pairs[mp_idx,1] == names(cresult))]]
         y              <- cresult[[which(model_pairs[mp_idx,2] == names(cresult))]]
         if( !is.na(attr(x,"model")$loglik) && !is.na(attr(y,"model")$loglik))
            cross_comparison_results[[cnames[mp_idx]]] <- cross_compare_models(x,y)
      }
   }
   return(cross_comparison_results)
}

`compute_pattern` <-
function(data,fun=mean){
   data                 <- data.matrix(aggregate(data,fun,by=list(class=data[,"class"])))[,-1]
   row.names(data)      <- 1:nrow(data)	
   if(!is.na(match("class", dimnames(data)[[2]])))
      data              <- data[,-match("class", dimnames(data)[[2]])]
   return(data)
}

`cross_compare_models` <-
function(x,y){
   # CONTINGENCY TABLE BY CROSS COMPARISONS BETWEEN MODEL 1 AND 2
   x_comparison	        <- as.table(ftable(attr(x,"model")$labelling,attr(y,"model")$labelling))
   # USE JOINT PROBABILITIES TO COMPUTE THE CONTINGENCY TABLE (%) 
   x_comparison_p       <- fullproba_ftable(attr(x,"model")$z,attr(y,"model")$z)
   map_rand_err_test    <- t.test(as.numeric(x_comparison/sum(x_comparison)-x_comparison_p))
   # CHI2_STATS: TEST THE ASSOCIATION BETWEEN THE TWO CATEGORICAL VARIABLES
   # H0: NO ASSOCIATION BETWEEN THE TWO VARIABLES
   # H1: THERE IS ASSOCIATION BETWEEN THE TWO VARIABLES
   ind_test             <- chisq.test(x_comparison)
   # CRAMER'S V 
   x_summary            <- summary(x_comparison)
   n                    <- as.numeric(x_summary["n.cases"]) 
   X2                   <- as.numeric(x_summary["statistic"]) 
   k                    <- min(dim(x_comparison)) 
   V                    <- sqrt(X2/(n*(k-1))) 
   # PEARSON'S CORRELATION: REPORT CORRELATION BETWEEN THE TWO CATEGORICAL
   # VARIABLES TEST (P) WHETHER THE CORRELATION IS EQUAL TO ZERO AND REPORTS
   # THE 95% CONFIDENCE INTERVAL
   cor_test     <- cor.test(attr(x,"model")$labelling,attr(y,"model")$labelling,use="pairwise.complete.obs")
   # BIND STATISTICS SUCH AS LOG-ODD RATIO, LAMBDA-SIBS INTO A FINAL TABLE 1RST
   # DIRECTION
   for(s in names(attr(x,"stats")))
      x_comparison	<- cbind(x_comparison,data.matrix(attr(x,"stats")[[s]][["out"]]))
   mat_stats	<- matrix(0,0,attr(y,"model")$G)
   for(s in names(attr(y,"stats")))
      mat_stats <- rbind(mat_stats, t(attr(y,"stats")[[s]][["out"]]))
   mat_stats    <- cbind(mat_stats, matrix(0,nrow(mat_stats),ncol(x_comparison)-attr(y,"model")$G))
   dimnames(mat_stats)[[2]] <- dimnames(x_comparison)[[2]]
   x_comparison <- rbind(x_comparison,mat_stats)
   # FINALLY, ORDER THE ROWS AND THE COLUMNS BY (DI)SIMILARIY, HIERARCHICAL
   # CLUSTERING
   ordering     <- list(row=hclust(dist(x_comparison[1:attr(x,"model")$G,1:attr(y,"model")$G]))$order,
                     column=hclust(dist(t(x_comparison[1:attr(x,"model")$G,1:attr(y,"model")$G])))$order)
   x_comparison	<- x_comparison[c(ordering$row,(attr(x,"model")$G+1):nrow(x_comparison)),c(ordering$col,(attr(y,"model")$G+1):ncol(x_comparison))]
   # BIND CHI2 STATS
   x_comparison <- rbind(x_comparison,"")
   x_comparison[which(x_comparison == 0)] <- ""
   row.names(x_comparison)[nrow(x_comparison)]    <- ind_test$method
   x_comparison[nrow(x_comparison),1:2]           <- c(ind_test$p.value,paste2(
                                                        "(X2=",ind_test$statistic,
                                                        ", df=",ind_test$parameter,")"))
   # BIND CRAMER'S V ASSOCIATION STAT
   x_comparison  <- rbind(x_comparison,"")
   row.names(x_comparison)[nrow(x_comparison)]    <- "Cramer's V"
   x_comparison[nrow(x_comparison),1]             <- V 
   # T-TEST TO EVALUATE WHETHER THE DIFFERENCES BETWEEN THE MAP CONTINGENCY
   # TABLE AND THE JOINT PROBABILITY DISTRIBUTION OF THE TWO FACTORS IS A
   # RANDOM ERROR OF MEAN = 0
   x_comparison  <- rbind(x_comparison,"")
   row.names(x_comparison)[nrow(x_comparison)]    <- "T-test: mapping error = 0?"
   x_comparison[nrow(x_comparison),1]             <- map_rand_err_test$p.value
   # TO QUICKEN THE LATER OPENING OF THE MANY CSV FILES,  WE REDUCE THE NUMBER
   # OF ZEROS TO CLEAR MANUALLY. IN WRITE.CSV2 'NA' ARE SUBSTITUTED BY EMPTY
   # STRINGS
   x_comparison[is.na(x_comparison) ]   <- ""
   x_comparison                         <- as.data.frame(x_comparison)
   return(x_comparison)
}

`fun_mbc_em` <-
function(data
        , G = 3
        , modelName = "EII"
        , rseed=6013) {
   # 
   G            <- as.numeric(G)
   modelName   <- as.character(modelName)
   # PREPARE A RANDOM MATRIX OF CLUSTER MEMBERSHIPS WITH G COLS AND N ROWS
   if(G >= 2){
      cmembership_matrix        <- matrix(rnorm(G*nrow(data)),nrow(data),G)
      # !!! NOT VERY SATISFYING THAT AT LEAST ONE OF THE COMPONENT IS ZERO...!!!
      cmembership_matrix        <- cmembership_matrix-apply(cmembership_matrix,1,min)
      cmembership_matrix        <- cmembership_matrix/apply(cmembership_matrix,1,sum)
   }
   else
      cmembership_matrix        <- matrix(1,nrow(data),G)
   # ESTIMATE AN INITIAL MODEL WHICH WE PROVIDE AS STARTING POINT TO EM
   m_step_estimate              <- mstep( data = data
                                        , modelName = modelName
                                        , z = cmembership_matrix)
   # AND THEN, STARTS EM WITH THE E-STEP
   model                       <- em(    data = data
                                        , modelName = m_step_estimate$modelName
                                        , parameters = m_step_estimate$parameters)
   #
   model[["labelling"]]        <- map(model[["z"]],warn=FALSE)
   # CALCULATE THE BIC SCORE
   model[["BIC"]]              <- bic(   modelName = modelName
                                        , loglik = model$loglik
                                        , n = nrow(data)
                                        , d = ncol(data)
                                        , G = G)
   return(model)
}

`stats_generalization` <-
function(data,class
   , fun_classifier = list(model=model_naive_bayes,predict=predict_naive_bayes)
   , K=7
   , fun_eval=stratified_traintest_split
   #, fun_eval=stratified_cv
   , title = NULL
)
{
   class                <- map(class)
   ldata                <- cbind(get_cdata(data),class=class)
   # INITIALIZATION
   indexes <- testfold  <- list()
   models <- predictions<- contingency.tables <- perfmeasures <- list()
   G                    <- length(table(class))
   if(!(G == 0)){
      # GENERATE THE DIFFERENT STRATA FROM THE K-FOLD CROSS VALIDATION
      out_cv <- fun_eval(ldata,K)
      indexes <- out_cv[["indexes"]]
      K <- out_cv[["K"]]
      # PROCEED TO EVALUATION FOR EACH K-FOLD
      for(k in 1:K){
         # INIT THE K-TH TEST FOLD
         testfold[[k]]          <- list()
         for(g in 1:G)
            testfold[[k]]       <- append(testfold[[k]],row.names(ldata[ldata[,"class"] == g,][indexes[[g]][k,],]))
         testfold[[k]]          <- unlist(testfold[[k]])
         # DEFINE TRAIN AND TEST SETS
         cdata_train            <- set_cdata(data = ldata[-pmatch(testfold[[k]], row.names(ldata)),]
                                                , settings = attr(data,"settings"))
         cdata_test             <- set_cdata(data = ldata[ pmatch(testfold[[k]], row.names(ldata)),]
                                                , tdata = attr(cdata_train,"tdata")
                                                , settings = attr(data,"settings"))
         trainset               <- cbind(get_cdata(cdata_train),class=cdata_train[,"class"])
         testset                <- cbind(get_cdata(cdata_test),class=cdata_test[,"class"])
         # TRAIN DIFFERENT MODELS, NAMELY NAIVE BAYES, SVM AND 1-NN RK: SVM
         # USES G(G-1) CLASSIFIERS, NAIVE BAYES AND KNN ARE SINGLE MULTI-CLASS
         models[[k]]            <- fun_classifier[["model"]](as.formula("class ~ ."),trainset,testset)
         predictions[[k]]       <- fun_classifier[["predict"]](models[[k]],testset)
         contingency.tables[[k]]<- table(predictions[[k]],testset[,"class"])
         perfmeasures[[k]]      <- sum(diag(contingency.tables[[k]]))/sum(contingency.tables[[k]])
      }
      # SUMMARY STATISTICS
      perfs		        <- unlist(perfmeasures)
      out                       <-t( matrix( c( mean(perfs),
                                                1.96*sd(perfs,na.rm=TRUE),K)
                                             , 3
                                             , G
                                             , dimnames=list(list(title,"CI95%","K-folds"),as.list(1:G)),byrow=FALSE))
      out[2:G,]                 <- NA
      return(list(      models = models
                        , predictions = predictions
                        , contingency.tables = contingency.tables
                        , perfmeasures = perfmeasures
                        , avg=mean(perfs,na.rm=TRUE)
                        , sd=1.96*sd(perfs,na.rm=TRUE)
                        , out=out))
   }
   else return(NULL)
}

`fullproba_ftable` <-
function(z1, z2){
   m_out                <- matrix(NA,ncol(z1),ncol(z2))
   row_name_vector      <- unique(c(row.names(z1),row.names(z2)))
   m                    <- matrix(NA, length(row_name_vector), ncol(z1)+ncol(z2),dimnames=list(row_name_vector,list()))
   m[row.names(z1),1:ncol(z1)] <- z1
   m[row.names(z2),(ncol(z1)+1):(ncol(z1)+ncol(z2))] <- z2
   for(i in 1:ncol(z1))
      for(j in 1:ncol(z2))
         m_out[i,j] <- mean(m[,i] * m[,j+ncol(z1)])
   return(m_out)
}

`generate_cdata_settings` <-
function(data){
   # DEFINE THE DEFAULT VALUES OF OUR CONF MATRIX
   conf_col_names <- list("group","in_canalysis","fun_transform","visu_groups","visu_ycoord")
   default_values <- c("",TRUE,"transform_AVG transform_SIGMA","var_group_1","NULL")
   # in_canalysis? (TRUE/FALSE) 
   #    Is the variable passed to the clustering algorithm?
   # visu_groups (CHAR VECTOR) 
   #    Group names, in use for the characterization by logodds and the
   #    visualization.
   # visu_ycoord? (INTEGER IN [1,10] OR NULL) 
   #    Shall we involve this variable in the parallel coordinates, if this is
   #    not set as null, then we report its values at the integer specified.
   conf_mat <- matrix(default_values,ncol(data),length(default_values),byrow=TRUE,
                  dimnames=list(colnames(data),conf_col_names))
   return(conf_mat)
}

`get_cdata` <-
function(data, which_data = NULL)
{
   var_names <- get_canalysis_variables(attr(data,"settings"))
   if(is.null(which_data))
      data_out <- data[,var_names]
   else
      data_out <- attr(data,"data_o")[,var_names]
   return(data_out)
}

`get_canalysis_variables` <-
function(settings)
{
   which_var         <- which(settings[,"in_canalysis"] == "TRUE")
   return(row.names(settings)[which_var])
}

`get_fun_stats` <-
function(fun_name="oddratios",...){
   # ODD RATIOS CHARACTERIZING EACH CLUSTER
   if(fun_name == "lambdasibs")
      return(function(data,class) { return(stats_lambdasibs(data,class,...)) })
   # ODD RATIOS CHARACTERIZING EACH CLUSTER
   if(fun_name == "oddratios")
      return(function(data,class) { return(stats_logodds(data,class,...)) })
   # AREA UNDER THE UNCERTAINTY CURVE
   if(fun_name == "auuc")
      return(function(data,class) { return(stats_auuc(class)) })
   # GENERALIZATION ESTIMATES OF MACHINE LEARNING ALGORITHMS
   if(fun_name == "gen_naive_bayes")
      return(function(data,class) { return(
               stats_generalization(
                  data,
                  class,
                  fun_classifier=list(model = model_naive_bayes,predict = predict_naive_bayes),
                  title = "naive Bayes",
                  ...)) })
   if(fun_name == "gen_knn")
      return(function(data,class) { return(
               stats_generalization(
                  data,
                  class,
                  fun_classifier=list(model = model_knn,predict = predict_knn),
                  title = "1 nearest neighbour",
                  ...)) })
   if(fun_name == "gen_svm")
      return(function(data,class) { return(
               stats_generalization(
                  data,
                  class,
                  fun_classifier=list(model = model_svm,predict = predict_svm),
                  title = "Support Vector Machine",
                  ...)) })
   else
      return(NULL)
}

`get_best_models` <-
function(cresult){
   bic_table    <- get_bic_table(cresult,which_bic = "relative_BIC", type="long")
   top_count    <- attr(cresult,"nbr_top_models")
   # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
   best_bics    <- bic_table[na.omit(order(bic_table$relative_BIC,decreasing=FALSE)[1:top_count]),]
   best_models  <- apply(best_bics[,c("modelName","G")],1,paste,collapse=",",sep="")
   return(best_models)
}

`init_data_cc` <-
function(data,settings)
{
   var_names         <- get_canalysis_variables(settings)
   data_mat          <- data[row.names(na.omit(data[,var_names])),]
   return(data_mat)
}


`get_bic_table` <- function(cresult,which_bic = "relative_BIC", type="wide"){
   bic_table    <- cbind(attr(cresult, "cfun_params"),BIC=NA,relative_BIC=NA)
   # PREPARE A BIC SCORE TABLE FROM THE RESULT
   for(i in 1:length(cresult))
      bic_table[i,"BIC"]                <- attr(cresult[[i]],"model")[["BIC"]]
   #
   bic_table[,"relative_BIC"]   <- 100*(bic_table[,"BIC"]/max(bic_table[,"BIC"],na.rm=TRUE)-1)
   #
   if(type == "wide"){
      bic_table                 <- aggregate(bic_table[,which_bic],by=list(modelName=bic_table$modelName,G=bic_table$G),FUN=mean)
      colnames(bic_table)       <- c("modelName","G",which_bic)
      bic_table                 <- reshape(bic_table[,c("G","modelName",which_bic)],idvar=c("G"),timevar=c("modelName"),direction="wide")
      row.names(bic_table)      <- bic_table[,"G"]
      colnames(bic_table)       <- gsub(paste2(which_bic,"."),"",colnames(bic_table))
      bic_table                 <- bic_table[,-1]
   }
   else
      bic_table               <- bic_table[,c("G","modelName",which_bic)]
   return(bic_table)
}



`plot.cdata` <- 
function(x, doPSOutput = TRUE, ...) {
   ps_output    <- paste2(attr(x,"prefix"),"_CANALYSIS_DATA_PLOT.ps")
   if(doPSOutput)       
      postscript(file=ps_output)
   ldata        <- get_cdata(x)
   ldata        <- as.data.frame(ldata[1:nrow(ldata),])
   par(mfrow = c(1,1) , new = FALSE)
   boxplot(ldata,'las'=2,main=paste2(ps_output," Boxplot "))
   par(mfrow = c(4,5))
   for(var in colnames(x)){
      # prepare tdata text
      var_idx   <- which(names(attr(x,"tdata")) == var)
      txt       <- paste2("(", length(table(x[,var])), " values\n")
      if(length(var_idx) == 1){
         for(tn in names(attr(x,"tdata")[[var_idx]])){
            if(tn == "transform_adjust")
               txt <-  paste2(txt, tolower(gsub("transform_","",tn)), "=", attr(x,"tdata")[[var_idx]][[tn]][["print"]], " ")
            else
               txt <-  paste2(txt, tolower(gsub("transform_","",tn)), "=", sprintf("%.2f",attr(x,"tdata")[[var_idx]][[tn]])," ")
            }
         }
      txt    <- paste2(txt,")")
      # plot the histogram
      hist(      x[,var]
               , main   = var
               , xlab   = txt 
               , ylab   = ""
               , cex    = 0.7)
   }
   if(doPSOutput)       
      graphics.off() 
}

`plot.cresult` <- 
function(x, doPSOutput = TRUE, query = NULL, ...) {
   if(is.null(query))   
      q_cmodel        <- 1:length(x)
   else 
      q_cmodel        <- query
   if(doPSOutput)
      postscript(file=paste2(attr(x,"prefix"),"_CRESULT_PLOT.ps"))
   for(i in q_cmodel){
      if(!is.na(attr(x[[i]],"model")$loglik)){
         par(mfrow=c(3,3), 'las'=1,'mai'=c(0.2,0.7,0.5,0.7),new=FALSE)
         for(fun_plot in unlist(attr(x,"fun_plot")))
            fun_plot(x[[i]])
      }
   }
   if(doPSOutput) 
      graphics.off() 
}

`print.cresult` <- 
function(x, html = TRUE, ...) {
   print_out <- list()
   print_out[["bic_table_relative"]]    <- get_bic_table(x,which_bic = "relative_BIC",type="wide")
   print_out[["bic_table_absolute"]]    <- get_bic_table(x,which_bic = "BIC",type="wide")
   print_out[["best_models"]]           <- get_best_models(x)
   print_xcomp                          <- compare_cresult(x)
   for(p_name in names(print_xcomp)){
      pxcomp                            <- as.matrix(print_xcomp[[p_name]])
      pxcomp2                           <- pxcomp
      pxcomp2[pxcomp2 == ""]            <- "0.00"
      for(i in 1:nrow(pxcomp)){
         for(j in 1:ncol(pxcomp)){
            if(!is.na(as.numeric(pxcomp[i,j]))){
               pxcomp2[i,j]             <- paste2("<div align='right'>",sprintf(pxcomp2[i,j],fmt='%.2f'),"</div>")
               if(as.numeric(pxcomp[i,j]) >= 1)
                  pxcomp2[i,j]             <- paste2("<font color='red'>",pxcomp2[i,j],"</font>")
               if(as.numeric(pxcomp[i,j]) <= -1)
                  pxcomp2[i,j]             <- paste2("<font color='blue'>",pxcomp2[i,j],"</font>")
            }
            else
               pxcomp2[i,j]             <- ""
         }
      }
      print_out[[p_name]]               <- structure(pxcomp2,csv=pxcomp)
   }
   if(html){
      directory <- getwd()
      file      <- paste2(attr(x,"prefix"),"_report.html")
      file      <- file.path(directory, file)
      cat("\n", file = file, append = TRUE)
      for(p_idx in 1:length(print_out)){
         HTML.title(names(print_out)[[p_idx]],file = file)
         file_csv       <- paste2(attr(x,"prefix"),"_",names(print_out)[p_idx],".csv")
         if(!is.null(attr(print_out[[p_idx]],"csv")))
            write.table(attr(print_out[[p_idx]],"csv"), na="", dec=",", sep=";", file=file_csv)
         else
            write.table(print_out[[p_idx]], na="", dec=",", sep=";", file=file_csv)
         HTML(paste2("<a href='",getwd(),"/",file_csv,"'>CSV file</a>"), file=file)
         HTML(print_out[[p_idx]],file=file,align="left")
      }
      for(m_name in names(x)){
         file_img       <- paste2(attr(x,"prefix"),"_",m_name,".png")
         img_size       <- c(1024,1024)
         png(file.path(directory,file_img),width=1.1*img_size[1],height=1.1*img_size[2],bg="transparent")
         plot(x,query = m_name,doPSOutput =FALSE)
         dev.off()
         HTMLInsertGraph(file_img,Caption = paste2("Figure. Graphic characteristics of model ",m_name)
               , file=file,GraphBorder = 0, WidthHTML=img_size[1],HeightHTML=img_size[2])
      }
      cat(paste2("Report written: ", file))
   }
   else
      return(print_out)
}

`fun_cmodel` <-
function(cresult) {
   # RETRIEVE DATA FOR THE CANALYSIS
   cdata        <- get_cdata(attr(cresult,"cdata"))
   for(i in 1:length(cresult)){
      cat("\n",names(cresult)[[i]])
      attr(cresult[[i]],"model")                <- cresult[[i]](cresult[[i]],cdata)
      if(!is.na(attr(cresult[[i]],"model")$loglik)){
         cdata_l                                <- cbind(attr(cresult,"cdata"),class=attr(cresult[[i]],"model")[["labelling"]])
         cat("-> Patterns")
         # COMPUTE PATTERNS
         for(p in names(attr(cresult,"fun_pattern")))
            attr(cresult[[i]],"pattern")[[p]]   <- compute_pattern(cdata_l,attr(cresult,"fun_pattern")[[p]])
         # FOR ORDERING, MAKE A DENDROGRAM FROM THE CENTER PATTERN (NUMBER 1)
         cat("-> Dendros")
         avg_p                                  <- attr(cresult[[i]],"pattern")[[1]]
         attr(cresult[[i]],"dendro_cluster")    <- hclust(dist(avg_p))
         attr(cresult[[i]],"dendro_var")        <- hclust(dist(t(avg_p)))
         # REORDER THE PATTERNS
         cat("-> Ordering")
         for(p in names(attr(cresult,"fun_pattern")))
            attr(cresult[[i]],"pattern")[[p]]   <- attr(cresult[[i]],"pattern")[[p]][
                                                                  attr(cresult[[i]],"dendro_cluster")$order
                                                                , rev(sort(colnames(avg_p)))]
         attr(cresult[[i]],"pattern")[["class_count"]] <-  table(cdata_l[,"class"])[attr(cresult[[i]],"dendro_cluster")$order]
         # SELECT DIVERGING COLORS FROM RColorBrewer
         attr(cresult[[i]],"cluster_colors")    <- brewer.pal(nrow(avg_p),"Set1")
         # STATS
         fun_stats                              <- attr(cresult,"fun_stats")
         cat("-> Stats")
         for(n in names(fun_stats))
            attr(cresult[[i]],"stats")[[n]]     <- fun_stats[[n]](data=attr(cresult,"cdata")
                                                        , class=attr(cresult[[i]],"model")[["z"]])
      }
   }
   # WRITE A BIC SCORE TABLE AS AN OUTPUT
   bic_table    <- get_bic_table(cresult,which_bic = "relative_BIC",type="wide")
   write.table(file=paste2(attr(cresult,"prefix"),"_Relative_BIC_Table.csv"),bic_table,dec=",",sep=";")
   return(cresult)
}

`model_knn` <-
function(formula, trainset, testset, k = 1, l = 0, prob = FALSE, use.all = TRUE){
   return(knn(    trainset[,-which(colnames(trainset) == "class")]
                , testset[ ,-which(colnames(testset) == "class")] 
                , trainset[,"class"]
                , k = k
                , l = l
                , prob = prob
                , use.all = use.all))
}

`model_naive_bayes` <-
function(formula, trainset, testset, laplace = 1){
   return(naiveBayes(formula, data=trainset, laplace = 1))
}

`model_svm` <-
function(formula, trainset, testset, kernel="linear",type="C-classification"){
   return(svm(formula, data=trainset,kernel=kernel,type=type))
}

`paste2` <-
function(...){
	return(paste(...,collapse="",sep=""))
}

`predict_knn` <-
function(model, testset){
   return(model)
}

`predict_naive_bayes` <-
function(model,testset){
   return(map(predict(model,testset[,-match("class",colnames(testset))],type="raw")))
}

`predict_svm` <-
function(model, testset){
   return(predict(model,testset[,-match("class",colnames(testset))]))
}

`write.cresult` <-
function(x, query = NULL){
   if(is.null(query))
      query             <- get_best_models(x)
   for(q in query){
      model             <- attr(x[[q]], "model")
      out               <- cbind(model[["z"]], class=map(model[["z"]]))
      colnames(out)[1:(ncol(out)-1)] <- 1:(ncol(out)-1)
      csv_output        <- paste2(attr(x,"prefix")
                                , "_Class_"
                                , model[["modelName"]]
                                ,"-"
                                ,model[["G"]]
                                ,".csv") 
      print(csv_output)
      write.csv2(out,file=csv_output)
   }
}

`set_cresult` <-
function( data          = NULL
        , cfun          = fun_mbc_em
        , cfun_settings = list(modelName=c("EII", "VII") ,G=1:9 ,rseed=6013)
        , fun_pattern   = list(mean=function(x){return(mean(x,na.rm=TRUE))}
                            , median=function(x){return(median(x,na.rm=TRUE))}
                            , lowquant=function(x){return(quantile(x,probs = 0.025,na.rm=TRUE))}
                            , upquant=function(x){return(quantile(x,probs = 0.975,na.rm=TRUE))})
        , fun_plot      = list(plot_parcoord=get_plot_fun(type="plot_parcoord")
                              , plot_legend=get_plot_fun(type="plot_legend")
                              , plot_image=get_plot_fun(type="plot_image")
                              , plot_dendro_cluster=get_plot_fun(type = "plot_dendro_cluster")
                              , plot_dendro_var=get_plot_fun(type="plot_dendro_var"))
        , fun_stats     = list(oddratios=get_fun_stats(fun_name="oddratios"))
        , nbr_top_models= 5
        ){
   # COMPUTES THE COMBINATION OF THE CLUSTERING FUNCTION PARAMETERS 
   param_set_df         <- cfun_settings[[1]]
   i <- 2
   while(i <= length(cfun_settings)){
      param_set_df      <- merge(param_set_df,cfun_settings[[i]])
      colnames(param_set_df)<-names(cfun_settings)[1:i]
      i <- i+1
   }

   # CONSTRUCT A LIST OF 'CMODEL' DATA STRUCTURES
   cresult_out          <- list()
   for(i in 1:nrow(param_set_df)){
      id_name           <- paste(as.matrix(param_set_df[i, names(cfun_settings)]),collapse=",")
      fun_call          <- function(cmodel,data){return(do.call(cfun,append(list(data=data),attr(cmodel,"cfun_settings"))))}
      cresult_out[[id_name]] <- structure(fun_call
                                , cfun_settings = param_set_df[i,]
                                , model         = NULL
                                , pattern       = list()
                                , dendro_cluster= NULL
                                , dendro_var    = NULL
                                , cluster_colors= NULL
                                , stats         = list()
                                , class         = "cmodel")
   }
   # CONSTRUCT A 'CRESULT' DATA STRUCTURE WITH CMODELS IN IT
   cresult_out          <- structure(cresult_out
                                , cdata   = data
                                , cfun_params      = param_set_df
                                , fun_plot         = fun_plot
                                , prefix           = attr(data,"prefix")
                                , fun_stats        = fun_stats
                                , fun_pattern      = fun_pattern
                                , nbr_top_models   = nbr_top_models
                                , class            = "cresult")
   return(cresult_out)
}

`set_cdata` <-
function(data, data_o = NULL, tdata = NULL, settings, init_fun = list(init_data_cc), prefix =""){
   tdata_list   <- list()
   settings     <- as.matrix(settings)
   # BACKUP DATA INTO DATA_O 
   if(is.null(data_o) )
      data_o <- data
   # DO INIT FUNCTIONS OF THE DATA
   for(fun in init_fun)
      data <- fun(data,settings)
   # TRANSFORM THE DATA
   t_result     <- transform_cdata(data,settings,tdata)
   data_out     <- data.matrix(t_result[["data"]])
   tdata_out    <- t_result[["tdata"]]
   cdata_out    <- structure(data_out
                        , data_o        = data_o  
                        , tdata         = tdata_out
                        , settings      = settings 
                        , prefix        = paste2(today(),"_",prefix)
                        , xlim          = c(-3,3)
                        , ylim          = c(0,50)
                        , class         = "cdata")
   return(cdata_out)
}

`plot_parcoord` <- 
function(xlim = c(-3,3)
      , cex = NULL
      , title = NULL
      , T_s = NULL
      , pattern = NULL){
   return(function(cmodel){
               plot(x=0 ,new=TRUE ,ann=FALSE ,pch=18,col="white",axes=FALSE ,
                     xlim=xlim , ylim=c(as.numeric(min(T_s[,"visu_ycoord"])),
                     as.numeric(max(T_s[,"visu_ycoord"]))))
               title(main = paste2("(",paste(as.matrix(attr(cmodel, "cfun_settings")),collapse=","),")",title), cex = cex * 0.8)
               for(x_name in pattern){
                  pattern         <- attr(cmodel,"pattern")[[x_name]][,row.names(T_s)]
                  lwd             <- 3
                  lty             <- "solid"
                  if(x_name != "mean") {
                     lwd          <- 1
                     lty          <- "dashed"
                  }
                  for(g in 1:nrow(pattern)){
                     D_i_s        <- cbind(X=as.numeric(pattern[g,]),Y=as.numeric(T_s[,"visu_ycoord"]))
                     D_i_s	 <- D_i_s[sort.list(D_i_s[,"Y"]),]
                     gap          <- D_i_s[2,"Y"] - D_i_s[1,"Y"]
                     for(l in 1:(nrow(D_i_s)-1))
                        if(!(D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap))
                           arrows(D_i_s[l,1], D_i_s[l,2], D_i_s[l+1,1], D_i_s[l+1,2], col=attr(cmodel,"cluster_colors")[g], length=0, lwd=lwd, lty=lty)
                        # ELSE, AS (is_white_gap == TRUE) THEN DO NOT DRAW ANY ARROW...
                  }
               }
               axis(2,at=as.numeric(T_s[,"visu_ycoord"]),labels=colnames(pattern),las=2,tick=FALSE)
               axis(1,at=seq(from = xlim[1] , to = xlim[2], by = (xlim[2] - xlim[1])/4))
         })
} 

`get_plot_fun` <-
function( type  = "plot_parcoord" 
        , cdata = NULL
        , title = NULL
        , xlim  = c(-3,3)
        , xy    = c(-2.2,0) 
        , cex   = 0.9
        , pattern ="mean"
        , color_gradient = brewer.pal(11,"RdBu")){
   if(type == "plot_parcoord"){
      cfun_parcoord     <- c()
      for(v in unique(as.character(attr(cdata,"settings")[,"visu_groups"]))){
         T_s            <- attr(cdata,"settings")[attr(cdata,"settings")[,"visu_groups"] == v,]
         lfun           <- eval(call("plot_parcoord",xlim = xlim , cex = cex , title = v, T_s = T_s , pattern = pattern))
         cfun_parcoord  <- c(cfun_parcoord, lfun)
      }
      return(cfun_parcoord)
   }
   if(type == "plot_image")            
      return(
         function(cmodel, apattern = pattern, acolgrad = color_gradient, atitle = title, acex = cex){
               x <- attr(cmodel,"pattern")[[apattern]]
               image(1L:nrow(x), 1L:ncol(x), x, xlim = 0.5 + c(0, nrow(x)), ylim = 0.5 + c(0, ncol(x)), axes = FALSE, 
                        xlab = "", ylab = "", col=acolgrad, main=atitle,add=FALSE)
               axis(2, 1L:ncol(x), labels = colnames(x), las = 2, line = -0.5, tick = 0, cex.axis = acex)
               axis(1, 1L:nrow(x), labels = row.names(x),las = 1, line = -0.5, tick = 0, cex.axis = acex)
            })
   if(type == "plot_legend")           
      return(
         function(cmodel, axy = xy){
            pattern              <- attr(cmodel, "pattern")[["class_count"]]
            for(g in 1:nrow(pattern)){
               g_name      <- attr(cmodel,"cluster_colors")[g]
               y1 <- y0    <- axy[2] + 5 * g / nrow(pattern)
               x0          <- axy[1]
               x1          <- 1.2 * axy[1]
               arrows(x0, y0, x1, y1, col=attr(cmodel,"cluster_colors")[g], length=0, lwd=3)
               text(x0+0.15, y1, labels=paste2(names(pattern)[g], " (",pattern[g],")"), pos=4)
            }
         })
   if(type == "plot_dendro_cluster" | type == "plot_dendro_var")   
        return( 
           function(cmodel
                ,  acex = cex
                , atitle = title ){
                     if(type == "plot_dendro_var") x <- attr(cmodel,"dendro_var")
                     else                            x <- attr(cmodel,"dendro_cluster")
                     plot(x, axes = FALSE, yaxs = "i", main=atitle, ylab=NULL,cex=acex)
           })
}

`set_tdata` <-
function(testimate=NULL,var=NULL){
  tdata_out <- structure(
       testimate
     , var   = var
     , class = "tdata")
  return(tdata_out)
}

`sibships_by_pairs` <-
function(data){
	# REMOVE INDIVIDUALS IN SIBSHIPS >= 3
	individual_list <- dimnames(data[data[,"family"]  %in% names(table(data[,"family"])[table(data[,"family"])>=2]) ,])[[1]]
	# FOR EACH FAMILY (ROW), LIST THE ASCERTAINED SIBLINGS (COLUMNS: 1,2,3,4,5)
	sibling_list		<- reshape(data[individual_list,c("UNIEK","family","member")],
							timevar="member", idvar="family", direction="wide")
	siblingpair_list	<- c()
	# FOR EACH FAMILY, CONSIDER ONLY THE TWO FIRST SIBLINGS
	for(s in 1:nrow(sibling_list))
		siblingpair_list<- rbind(siblingpair_list,t(sibling_list[s,!is.na(sibling_list[s,])][2:3]))
	return(list(data=data[siblingpair_list,],model=NULL,tfun="sibpairs"))
}

`stats_lambdasibs` <-
function(data,class)
{
   class <- map(class)
   b                    <- cbind(data[["orig"]],class = class)
   b                    <- reshape(b[,c("family","member","class")],timevar="member",idvar="family",direction="wide")
   b                    <- b[,-1]
   # MAKE A DATA STRUCTURE FROM PREVIOUSLY RESHAPED 'B' SPARSE MATRIX 
   # . IN THE 1ST COLUMN THE PROBANT, 
   # . IN THE 2ND ITS SIBLING 
   # . ALL THE NA'S ARE REMOVED
   for(j in 1:2)
      for(i in 1:nrow(b))
         while(is.na(b[i,j]))
            b[i,j:(ncol(b)-1)] <- b[i,(j+1):ncol(b)]
   b <- b[,1:2]
   class_set            <- levels(factor(as.matrix(class)))
   sibship_size		<- dim(b)[[1]]
   rowNames		<- c("Counts","P_cl","var_cl","P_cl_s1","P_cl_s2","P_cl_s1_s2",
                                          "var_cl_s1_s2","cov_P_cl_P_cl_s1_s2","lambdaSib","var_lambda",
                                          "se_lambda","ci_lambda","Significant?")
   s                    <- as.data.frame(matrix(0,length(rowNames),length(class_set)+1))
   # SET DIM NAMES FOR THE RESULT TABLE 'S'
   dimnames(s)[[1]]     <- rowNames
   dimnames(s)[[2]]     <- c(class_set,"Total")
   # COUNT THE NUMBER OF INDIVIDUALS IN EACH CLASS
   s["Counts",names(table(class))] <- table(class)
   s["Counts","Total"]  <- sum(s["Counts",names(table(class))])
   # P_cl          := CLASS PREVALENCE
   s["P_cl",]           <- s["Counts",]/s["Counts","Total"]
   # var_cl        := VAR( CLASS PREVALENCE ESTIMATE)
   s["var_cl",]         <- s["P_cl",]*(1-s["P_cl",])/s["Counts","Total"]
   # P_cl_s1       := CLASS PREVALENCE FOR SIB_1
   s["P_cl_s1",]        <- c(table(factor(b[,1],levels=levels(as.factor(class)))),sibship_size)/sibship_size
   # P_cl_s2       := CLASS PREVALENCE FOR SIB_2
   s["P_cl_s2",]        <- c(table(factor(b[,2],levels=levels(as.factor(class)))),sibship_size)/sibship_size
   # P_cl_s1_s2    := P ( s_i^1 AND s_i^2 )
   for(i in as.numeric(class_set))
            s["P_cl_s1_s2",as.character(i)] <- dim(b[b[,1] == i & b[,2] == i,])[[1]]/sibship_size
   # var_cl_s1_s2  := VAR (P ( s_i^1 AND s_i^2 ))
   s["var_cl_s1_s2",] <- s["P_cl_s1_s2",]*(1-s["P_cl_s1_s2",])/sibship_size
   # lambdasib     := ESTIMATION OF THE RISK RATIO
   s["lambdaSib",]    <- s["P_cl_s1_s2",]/s["P_cl",]^2
   # "cov_P_cl_P_cl_s1_s2",
   s["cov_P_cl_P_cl_s1_s2",] <- s["P_cl_s1_s2",]*(1-s["P_cl",])/s["Counts","Total"]
   # var_lambda
   s["var_lambda",]   <- (1/s["P_cl",]^4)*(s["var_cl_s1_s2",]
                                    -4*s["cov_P_cl_P_cl_s1_s2",]*s["P_cl_s1_s2",]/s["P_cl",]
                                    +4*s["var_cl",]*(s["P_cl_s1_s2",]/s["P_cl",])^2)
   # se_lambda
   s["se_lambda",]   <- sqrt(s["var_lambda",]) 
   # ci_lambda
   s["ci_lambda",]   <- 1.96*s["se_lambda",]
   s["Significant?",]<- (s["lambdaSib",]-s["ci_lambda",])>1
   s["Significant?","Total"]<- sum(s["Significant?",1:length(class_set)])
   #
   return(list(out=t(s[,1:length(class_set)])))
}

`stats_logodds` <-
function(data,class,fun_midthreshold=median){
   #
   class  <- map(class)
   ldata  <- cbind(data,class=class)
   sgroup <- as.matrix(attr(data,"settings")[,"group"])
   s      <- as.data.frame(matrix(0,length(unique(class)),length(unique(sgroup[,1]))))
   # FOR EACH HIERARCHICAL SUBSET OF OUTCOMES, MAKE SUM SCORES
   for(gr in unique(sgroup[,1])){
      mat.l   <- cbind(SScore=apply(data.matrix(ldata[,row.names(sgroup)[sgroup[,1] == gr]]),1,sum),class=class)
      # FOR EACH PATTERN, DERIVE ITS STATISTIC
      for(g in sort(unique(class))){
         med_SScore <- fun_midthreshold(mat.l[,"SScore"])
         m11 <- nrow(mat.l[mat.l[,"class"]==g & mat.l[,"SScore"] >= med_SScore,])
         if(is.null(m11))
                  m11 <- 0
         m12 <- nrow(mat.l[mat.l[,"class"]==g & mat.l[,"SScore"] < med_SScore,])
         if(is.null(m12))
                  m12 <- 0			
         m21 <- nrow(mat.l[mat.l[,"class"]!=g & mat.l[,"SScore"] >= med_SScore,])
         if(is.null(m21))
                  m21 <- 0			
         m22 <- nrow(mat.l[mat.l[,"class"]!=g & mat.l[,"SScore"] < med_SScore,])
         if(is.null(m22))
                  m22 <- 0			
         mat <- matrix(c(m11,m12,m21,m22),2,2,byrow=TRUE)
         # LOG OF THE ODD-RATIO ALSO NAMED CROSS PRODUCT
         s[g,match(gr,unique(sgroup))] <- log(mat[1,1] * mat[2,2] / (mat[1,2] * mat[2,1]))
         dimnames(s)[[2]][match(gr,unique(sgroup))] <- paste(gr,sprintf("%.1f",med_SScore),sep="_",collapse="")
      }
   }
   return(list(out=s))
}


`stratified_traintest_split` <-
function(data,K=10,perc=0.7,rseed = 6013){
   indexes <- list()
   G <- length(table(data[,"class"]))
   for(g in 1:G){
      strata_idx <- nrow(data[data[,"class"] == g,])
      indexes[[g]] <- matrix(0,0,strata_idx * (1-perc)) 
      set.seed(rseed)
      for(k in 1:K)
            indexes[[g]] <- rbind(indexes[[g]],sample(1:strata_idx)[1:(strata_idx * (1-perc))])
   }
   return(list(indexes=indexes,K=K))
}

`today` <-
function(){
   tmp_date <- format(Sys.time(), "%Y-%m-%d")
   return(tmp_date)
}

`transform_ABSMAX` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(max(abs(x),na.rm=TRUE))}
                      , function(x){return(do.call('/',x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_ALL` <-
function(data,tfun=list(function(x){return(0)},function(x){return(do.call('+',x))}),tdata){
   if(is.null(attr(tdata,"tfun"))) 
      testimate         <- tfun[[1]](data)
   else 
      testimate         <- as.numeric(tdata)
   vdata                <- tfun[[2]](list(data,testimate))
   tdata_out            <- set_tdata(testimate=testimate,var=attr(tdata,"var"))
   return(list(data=vdata,tdata=tdata_out))
}


`transform_AVG` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(mean(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_L1` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun      <- list(function(x){return(x/sqrt(sum(x,na.rm=TRUE)))}
                      , function(x){return(do.call('/',x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_L2` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun      <- list(function(x){return(x/sqrt(sum(x^2,na.rm=TRUE)))}
                      , function(x){return(do.call('/',x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_MAX` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(max(x,na.rm=TRUE))}
                      , function(x){return(do.call('-', x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_MEDIAN` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(median(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_MIN` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun       <- list(function(x){return(min(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}

`transform_SIGMA` <-
function(data,tdata){
   vdata        <- data[,attr(tdata,"var")]
   tfun         <- list(function(x){return(sd(x,na.rm=TRUE))}
                      , function(x){return(do.call('/', x))})
   return(transform_ALL(data=vdata,tfun=tfun,tdata=tdata))
}


`transform_adjust` <-
function(data, tdata, tformula){
   var          <- strsplit(tformula,"~")[[1]][1]
   data         <- as.data.frame(data)
   if(length(tdata) == 0){
      print(tformula)
      tdata_lm          <- lm(as.formula(tformula),data=data)
      tdata_lm[["print"]]<- tformula
      tdata             <- set_tdata(tdata_lm,var=attr(tdata,'var'))
      # AS THE VARS WHICH ARE NOT INVOLVED IN THE CANALYSIS MAY HAVE MISSING
      # VALUES, TO PRESERVE THE NA'S WE REBUILD A MATRIX DATA STRUCTURE WHICH
      # SERVES TO INDEX PROPERLY THE RESIDUALS. AND WE RETURN A SIMPLE NUMERIC 
      # VECTOR WHICH HAS SOME NA'S IN IT.
      vdata             <- matrix(NA,nrow(data),1,dimnames=list(row.names(data),var))
      tresiduals        <- residuals(tdata)
      vdata[names(tresiduals),] <- tresiduals
      vdata             <- as.numeric(vdata)
   }
   else
      vdata             <- data[,var]-predict(tdata,data)
   return(list(data=vdata,tdata=tdata))
}

`transform_cdata` <-
function(data,settings,tdata=NULL){
   if(is.null(tdata)) tdata <- list()
   # FOR EACH VAR
   for(var in row.names(settings)){
      tdata_idx                 <- which(names(tdata) == var)
      if(length(tdata_idx) == 0) 
         tdata[[var]]  <- list()
      # USE (MULTIPLE)-SPACES AS SEPARATOR BETWEEN TRANSFORMATION FUNCTIONS
      fun_transform             <- strsplit(settings[var,"fun_transform"],"[ ]+")[[1]]
      for(lfun in fun_transform){
         # IN CASE 'LFUN' AS AN INSIDE PARAMETER, RETRIEVE IT (1) REPLACE ALL
         # ('\") BY SPACES AND THEN (2) SPLIT
         lfun                   <- strsplit(gsub("[('\")]"," ",lfun),"[ ]+")[[1]]
         if(length(tdata_idx) == 0)
            tdata_var_lfun      <- set_tdata(testimate=NULL,var=var)
         else
            tdata_var_lfun      <- tdata[[var]][[lfun[1]]]
         arg_list               <- list(data,tdata=tdata_var_lfun)
         # IF LFUN HAS SOME ARGUMENTS, WE APPEND THEM TO THE ARGUMENT LIST
         if(!is.na(lfun[2]))
            arg_list            <- append(arg_list,lfun[2])
         # PROCEED TO THE TRANSFORM AND OUTPUT A DATA VECTOR AND A TDATA
         var_result             <- do.call(lfun[1],args=arg_list)
         # UPDATE THE MAIN DATA FRAME WITH THE VECTOR
         data[,var]             <- var_result[["data"]]
         # KEEP RECORD OF EACH INDIVIDUAL TRANSFORMATION FOR THE 'VAR'
         tdata[[var]][[lfun[1]]] <- var_result[["tdata"]]
      }
   }
   return(list(data=data,tdata=tdata))
}

`transform_sibordering` <-
function(data){
   for(f in unique(data$family)){
      fmembers <- sort(data[data$family == f,"member"])
      if(!is.na(fmembers[1]) && !is.na(fmembers[2])){
         data[data$family == f & data$member == fmembers[1],"member"] <- 1
         data[data$family == f & data$member == fmembers[2],"member"] <- 2 
      }
   }
   return(data[data$member <= 2,])
}

