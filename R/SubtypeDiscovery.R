`analysis` <-
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
   cresult              <- cresult_numeric_characteristics(cresult, data, config)
   # COMPARE BEST CLUSTERS  2-BY-2 AND SAVE OUTPUT INTO CSV FILES
   compare_canalysis(cresult, config=config)
   return(cresult)
}

`area_under_uncertainty_curve` <-
function(class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}

`compare_canalysis` <-
function(analysis1,analysis2=NULL,query=NULL,config){
   prefix               <- attr(config,"prefix")
   bbox_threshold       <- attr(config,"bbox_threshold")
   cross_comparison_results <- model <- list()
#   if(!is.null(analysis2) && !is.null(query)){
#      # WE COMPARE THE RESULT OF TWO ANALYSIS FOR A GIVEN MODEL
#      canalysis2 <- analysis2[["cluster_analysis"]]
#      model[[1]] <- retrieve_model(canalysis1,query)
#      model[[2]] <- retrieve_model(canalysis2,query)
#      if(!is.null(filePrefix))
#         filePrefix <- paste2(filePrefix,"_X_comparison_",gsub(",","_",query),".csv")
#      if(!is.na(model[[1]][["model"]]$loglik) && !is.na(model[[2]][["model"]]$loglik))
#         cross_comparison_results[[1]] <- cross_compare_models(model[[1]],model[[2]],fileNameCSV=filePrefix)
#   }
#   else{
   # SELECT ALL 2-BY-2 COMPARISONS BETWEEN MODELS THAT ARE WITHIN THE RELATIVE
   # 5% BBOX TO THE BEST BIC VALUE
   vertices <- get_vertice_set(analysis1[["bic_table_relative"]], bbox_threshold)
   if(nrow(vertices)>=1){
      # PROCEED TO 2-BY-2 COMPARISON OVER THE SET OF VERTICES
      for(v in 1:nrow(vertices)){
         # RETRIEVE CLUSTERING Z-MAPS FROM 'CANALYSIS'
         model[[1]] <- retrieve_model(analysis1,vertices[v,1])
         model[[2]] <- retrieve_model(analysis1,vertices[v,2])
         lprefix <- paste2(prefix,"_",gsub(",","_",vertices[v,1]),"-",gsub(",","_",vertices[v,2]),".csv")
         if( !is.na(model[[1]][["model"]]$loglik) && !is.na(model[[2]][["model"]]$loglik))
            cross_comparison_results[[v]] <- cross_compare_models(model[[1]],model[[2]],fileNameCSV=lprefix)
      }
   }
   for(m in unique(as.character(vertices))){
      m <- retrieve_model(analysis1,m)
      save_model(m[["model"]],attr(config,"prefix"))
   }

#   }
   return(cross_comparison_results)
}

`compute_plot_data` <-
function(data,fun=mean){
   data                 <- data.matrix(aggregate(data,fun,by=list(class=data[,"class"]),na.rm=TRUE))[,-1]
   row.names(data)      <- 1:nrow(data)	
   if(!is.na(match("class", dimnames(data)[[2]])))
      data              <- data[,-match("class", dimnames(data)[[2]])]
   return(data)
}

`cresult_graphic_characteristics` <-
function(cresult,canalysis_data,canalysis_config){
   # RETRIEVE VIGNETTE TITLES 
   parcoord_titles      <- unique(as.character(attr(canalysis_data,"settings")[,"var_by_groups"]))
   fig_params           <- list()
   for(v in parcoord_titles)
      fig_params        <- c(fig_params, list(list(title = v , type = "plot_parcoord")))
   fig_params           <- c(fig_params, list(list(title = "", type = "plot_legend")))
   fig_params           <- c(fig_params, list(list(title = "", type = "plot_image")))
   fig_params           <- c(fig_params, list(list(title = "", type = "plot_dendro_cluster")))
   fig_params           <- c(fig_params, list(list(title = "", type = "plot_dendro_var")))
   # PAGE-MARGINS
   par(mfrow=c(3,3), 'las'=1,'mai'=c(0.4,0.55,0.18,0.05))
   for( m_idx in 1:length(cresult[["out"]]) ){
      m         <- cresult[["out"]][[m_idx]]
      pdata     <- list()
      if(!is.na(m[["loglik"]]) & m[["G"]] >= 3){
         for(fig in fig_params)
            pdata <- c(pdata, list(set_plot_data(canalysis_data = canalysis_data , 
                        canalysis_config = canalysis_config , model = m , 
                        plot_data_fun = list(center_pattern  = mean) , type = fig[["type"]], 
                        title_str = fig[["title"]])))
         # STORE THE PDATA INTO THE MODEL
         cresult[["out"]][[m_idx]][["pdata"]] <- pdata
         # DO THE PLOTTING
         for(local_pdata in pdata)
            attr(local_pdata,"plot_fun")(local_pdata)
         }
   }
   return(cresult)
}

`cresult_numeric_characteristics` <-
function(cresult,canalysis_data,canalysis_config){
   for( m_idx in 1:length(cresult[["out"]]) ){
      m         <- cresult[["out"]][[m_idx]]
      pdata     <- list()
      if(!is.na(m[["loglik"]]) & m[["G"]] >= 3){
         list_stat_fun <- attr(canalysis_config,"stats_fun")
         for(n in names(list_stat_fun)){
            cresult[["out"]][[m_idx]][["stats"]][[n]] <- list_stat_fun[[n]](data=canalysis_data,class=m[["z"]])
            cat(m[["G"]],m[["modelName"]],", ")
         }
      }
   }
   return(cresult)
}

`cross_compare_models` <-
function(model1,model2,fileNameCSV=today()){
   # CONTINGENCY TABLE BY CROSS COMPARISONS BETWEEN MODEL 1 AND 2
   x_comparison	  <- as.table(ftable(model1[["labelling"]],model2[["labelling"]]))
   # USE JOINT PROBABILITIES TO COMPUTE THE CONTINGENCY TABLE (%) 
   x_comparison_p <- fullproba_ftable(model1[["model"]]$z,model2[["model"]]$z)
   map_rand_err_test <- t.test(as.numeric(x_comparison/sum(x_comparison)-x_comparison_p))
   # CHI2_STATS: TEST THE ASSOCIATION BETWEEN THE TWO CATEGORICAL VARIABLES
   # H0: NO ASSOCIATION BETWEEN THE TWO VARIABLES
   # H1: THERE IS ASSOCIATION BETWEEN THE TWO VARIABLES
   ind_test     <- chisq.test(x_comparison)
   # CRAMER'S V 
   x_summary    <- summary(x_comparison)
   n            <- as.numeric(x_summary["n.cases"]) 
   X2           <- as.numeric(x_summary["statistic"]) 
   k            <- min(dim(x_comparison)) 
   V            <- sqrt(X2/(n*(k-1))) 
   # PEARSON'S CORRELATION: REPORT CORRELATION BETWEEN THE TWO CATEGORICAL
   # VARIABLES TEST (P) WHETHER THE CORRELATION IS EQUAL TO ZERO AND REPORTS
   # THE 95% CONFIDENCE INTERVAL
   cor_test     <- cor.test(model1[["labelling"]],model2[["labelling"]],use="pairwise.complete.obs")
   # BIND STATISTICS SUCH AS LOG-ODD RATIO, LAMBDA-SIBS INTO A FINAL TABLE 1RST
   # DIRECTION
   for(s in names(model1[["model"]][["stats"]]))
      x_comparison	<- cbind(x_comparison,data.matrix(model1[["model"]][["stats"]][[s]][["out"]]))
   mat_stats	<- matrix(0,0,model2[["G"]])
   for(s in names(model2[["model"]][["stats"]]))
      mat_stats <- rbind(mat_stats, t(model2[["model"]][["stats"]][[s]][["out"]]))
   mat_stats <- cbind(mat_stats, matrix(0,nrow(mat_stats),ncol(x_comparison)-model2[["G"]]))
   dimnames(mat_stats)[[2]] <- dimnames(x_comparison)[[2]]
   x_comparison <- rbind(x_comparison,mat_stats)
   # FINALLY, ORDER THE ROWS AND THE COLUMNS BY (DI)SIMILARIY, HIERARCHICAL
   # CLUSTERING
   ordering <- list(row=hclust(dist(x_comparison[1:model1$G,1:model2$G]))$order,
                                 column=hclust(dist(t(x_comparison[1:model1$G,1:model2$G])))$order)
   x_comparison	<- x_comparison[c(ordering$row,(model1$G+1):nrow(x_comparison)),c(ordering$col,(model2$G+1):ncol(x_comparison))]
   # BIND CHI2 STATS
   x_comparison  <- rbind(x_comparison,"")
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
   x_comparison[is.na(x_comparison) ] <- ""
   if(!is.null(fileNameCSV)){
      # WRITE OUTPUT IN A CSV FILE
      print(fileNameCSV)
      write.table(x_comparison,file=fileNameCSV,na="",dec=",",sep=";")
   }
   return(x_comparison)
}

`em_clustering` <-
function(data, G_set=1:9,
        modelNames=c("EII", "VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"
        ),rseed=6013) {
   set.seed(rseed)
   # clustering results
   cresult_tmp      <- list()
   idx          <- 1
   BIC          <- matrix(NA,length(G_set),length(modelNames),dimnames=list(G_set, modelNames))
   for(g in G_set){
      for(m in modelNames){
         cresult_tmp[[idx]]             <- em_mbc(data=data,modelName=m, G=g)
         BIC[as.character(g),m]         <- bic(modelName = m, loglik = cresult_tmp[[idx]]$loglik,n=nrow(data),d=ncol(data),G=g)
         idx                            <- idx+1
      }
   }
   best_model_idx       <- which(BIC == max(BIC,na.rm=TRUE))[1]
   cresult              <- cresult_tmp[[best_model_idx]]
   cresult[["out"]]     <- cresult_tmp
   attr(cresult[["out"]],"BIC") <- BIC
   return(cresult)
}

`em_mbc` <-
function(data, modelName="VVI", G=5) {
   if(G >= 2){
      mat          <- matrix(rnorm(G*nrow(data)),nrow(data),G)
      # not very satisfying that at least one of the component is zero...
      mat          <- mat-apply(mat,1,min)
      mat          <- mat/apply(mat,1,sum)
   }
   else
      mat          <- matrix(1,nrow(data),G)
   #
   msEst        <- mstep(modelName = modelName, data=data,z=mat)
   # an initial model is then estimated and provided as starting point to the em-algorithm
   em_model     <- em(modelName = msEst$modelName, data = data,parameters = msEst$parameters)
   return(em_model)
}

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

`generate_canalysis_config` <-
function(data){
   # DEFINE THE DEFAULT VALUES OF OUR CONF MATRIX
   conf_col_names <- list("group","in_canalysis","fun_transform","var_by_groups","visu_ycoord")
   default_values <- c(TRUE,"transform_AVG transform_SIGMA","var_group_1","NULL")
   # in_canalysis? (TRUE/FALSE) 
   #    Is the variable passed to the clustering algorithm?
   # var_by_groups (CHAR VECTOR) 
   #    Group names, in use for the characterization by logodds and the
   #    visualization.
   # visu_ycoord? (INTEGER IN [1,10] OR NULL) 
   #    Shall we involve this variable in the parallel coordinates, if this is
   #    not set as null, then we report its values at the integer specified.
   conf_mat <- matrix(default_values,ncol(data),length(default_values),byrow=TRUE,
                  dimnames=list(colnames(data),conf_col_names))
   return(conf_mat)
}

`get_canalysis_data` <-
function(data, which_data = NULL)
{
   if(class(data) != "canalysis_data")
   {
      print("Error (get_canalysis_data): You must provide data of type canalysis_data")
      break ;
   }
   else
   {
      var_names <- get_canalysis_variables(attr(data,"settings"))
      if(is.null(which_data))
         data_out <- data[,var_names]
      else
         data_out <- attr(data,"data_o")[,var_names]
      return(data_out)
   }
}

`get_canalysis_variables` <-
function(settings)
{
   which_var         <- which(settings[,"in_canalysis"] == "TRUE")
   return(row.names(settings)[which_var])
}

`get_coloring_scheme` <-
function(data){
   # COMPUTE THE PATTERN, 1ST
   data_pattern         <- compute_plot_data(as.matrix(data),mean)
   # MAKE A DENDROGRAM FROM THESE AVG PATTERN FOR THE DENDROGRAM ORDERING
   dendro_patterns      <- hclust(dist(data_pattern))
   # SELECT DIVERGING COLORS FROM RColorBrewer
   cluster_colors       <- brewer.pal(nrow(data_pattern),"Set1")
   # RETURN A DATA STRUCTURE WITH CLUSTERS RIGHTLY ORDERED AND THEIR ASSOCIATED COLOR
   return(list(cluster_order=dendro_patterns$order,cluster_color=cluster_colors))
}

`get_model_set` <-
function(relativeBic, bbox_threshold=5){
   relativeBic <- list(bicTable=relativeBic, bbox = apply(abs(relativeBic) < bbox_threshold,2,which))
   relativeBic$bbox <- relativeBic$bbox[which(lapply(relativeBic$bbox,length) >0)]
   # DERIVE VERTICES TO COMPARE WITHIN THE 5% BBOX OF THE BIC TABLE
   relativeBic$best_models <- matrix(0,0,2,dimnames=list(list(),list("modelName","G")))
   for(m in names(relativeBic$bbox))
      for(g in relativeBic$bbox[[m]])
         relativeBic$best_models <- rbind(relativeBic$best_models,c(m,row.names(relativeBic[["bicTable"]])[g]))
   return(relativeBic[["best_models"]])
}

`get_stats_fun` <-
function(fun_name="oddratios",...){
   # ODD RATIOS CHARACTERIZING EACH CLUSTER
   if(fun_name == "oddratios")
      return(
         function(data,class) { return(statistics_logodds(data,class,...)) })
   # AREA UNDER THE UNCERTAINTY CURVE
   if(fun_name == "auuc")
      return(
         function(data,class) { return(area_under_uncertainty_curve(class)) })
   # GENERALIZATION ESTIMATES OF MACHINE LEARNING ALGORITHMS
   if(fun_name == "gen_naive_bayes")
      return(
         function(data,class) { return(
               evaluate_generalization(
                  data,
                  class,
                  fun_classifier=list(model = model_naive_bayes,predict = predict_naive_bayes),
                  title = "naive Bayes",
                  ...)) })
   if(fun_name == "gen_knn")
      return(
         function(data,class) { return(
               evaluate_generalization(
                  data,
                  class,
                  fun_classifier=list(model = model_knn,predict = predict_knn),
                  title = "1 nearest neighbour",
                  ...)) })
   if(fun_name == "gen_svm")
      return(
         function(data,class) { return(
               evaluate_generalization(
                  data,
                  class,
                  fun_classifier=list(model = model_svm,predict = predict_svm),
                  title = "Support Vector Machine",
                  ...)) }
         )
   else
      return(NULL)
}

`get_vertice_set` <-
function(relativeBic, bbox_threshold=5){
   relativeBic <- list( bicTable=relativeBic, 
                        # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
                        bbox = apply(relativeBic <= max(sort(relativeBic)[1:bbox_threshold],na.rm=TRUE), 2,which)
                        # TO SELECT ALL MODELS BELOW A GIVEN THRESHOLD
#                        bbox = apply(abs(relativeBic) < bbox_threshold,2,which)
                        )
   relativeBic$bbox <- relativeBic$bbox[which(lapply(relativeBic$bbox,length) >0)]
   # DERIVE VERTICES TO COMPARE WITHIN THE 5% BBOX OF THE BIC TABLE
   relativeBic$vertices <- matrix(0,0,2,dimnames=list(list(),list("m1_g1","m2_g2")))
   for(m1 in names(relativeBic$bbox))
      for(m2 in names(relativeBic$bbox))
         if(match(m2,names(relativeBic$bbox)) >= match(m1,names(relativeBic$bbox)))
            for(g1 in relativeBic$bbox[[m1]])
               for(g2 in relativeBic$bbox[[m2]])
                  if((m1 != m2 & g2 >= g1) | (m1 == m2 & g2 > g1))
                     relativeBic$vertices <- rbind(relativeBic$vertices,c(
                                                      paste(m1,row.names(relativeBic[["bicTable"]])[g1],collapse="",sep=","),
                                                      paste(m2,row.names(relativeBic[["bicTable"]])[g2],collapse="",sep=",")))
   return(relativeBic[["vertices"]])
}

`init_data_cc` <-
function(data,settings)
{
   var_names         <- get_canalysis_variables(settings)
   data_mat          <- data[row.names(na.omit(data[,var_names])),]
   return(data_mat)
}

`mclust2` <-
function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = NULL, warn = FALSE, ...) 
{
    mc                  <- match.call(expand.dots = FALSE)
    mc[[1]]             <- as.name("mclust2_BIC")
    Bic                 <- eval(mc, parent.frame())
    G                   <- attr(Bic, "G")
    modelNames          <- attr(Bic, "modelNames")
    Sumry               <- summary_mclust2_BIC(Bic, data, G = G, modelNames = modelNames)
    if (!(length(G) == 1)) {
        bestG           <- length(unique(Sumry$cl))
        if (bestG == max(G)) 
            warning("optimal number of clusters occurs at max choice")
        else if (bestG == min(G)) 
            warning("optimal number of clusters occurs at min choice")
    }
    attr(Bic, "n")      <- attr(Bic, "warn") <- NULL
    attr(Bic, "initialization") <- attr(Bic, "control") <- NULL
    attr(Bic, "d")      <- attr(Bic, "returnCodes") <- attr(Bic, "class") <- NULL
    oldClass(Sumry)     <- NULL
    Sumry$bic           <- Sumry$bic[1]
    ans                 <- c(list(out = Bic), Sumry)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty","out")
    structure(ans[orderedNames], class = "mclust2_mclust2_Mclust")
}

`mclust2_BIC` <-
function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), 
    initialization = list(hcPairs = NULL, subset = NULL, noise = NULL), 
    Vinv = NULL, warn = FALSE, x = NULL, ...) 
{
    #
    # CONTROL OF THE PARAMETERS
    # 
    if (!is.null(x)) {
        if (!missing(prior) || !missing(control) || !missing(initialization) || !missing(Vinv)) 
            stop("only G and modelNames may be specified as arguments when x is supplied")
        #
        # SET LOCAL SCOPE VARIABLES WITH ATTRIBUTES OF THE SAME NAME FROM OBJECT x
        #
        prior           <- attr(x, "prior")
        control         <- attr(x, "control")
        initialization  <- attr(x, "initialization")
        Vinv            <- attr(x, "Vinv")
        warn            <- attr(x, "warn")
    }
    #
    # DATA MATRIX CONTROL (DIMENSION)
    #
    dimData             <- dim(data)
    oneD                <- is.null(dimData) || length(dimData[dimData > 1]) == 1
    if (!oneD && length(dimData) != 2) 
        stop("data must be a vector or a matrix")
    # SET LOCAL SCOPE DIMENSION VARIABLES
    if (oneD) {
        data            <- drop(as.matrix(data))
        n               <- length(data)
        d               <- 1
    }
    else {
        data            <- as.matrix(data)
        n               <- nrow(data)
        d               <- ncol(data)
    }
    #
    # SET DEFAULT MODEL PARAMETERS IF NONE ARE PROVIDED (X IS NULL)
    #
    if (is.null(x)) {
        if (is.null(modelNames)) {
            # IF DIMENSION IS OF LENGTH 1, ONLY TWO AVAILABLE MODELS
            if (d == 1) {
                modelNames <- c("E", "V")
            }
            # MODELS FOR OTHER DIMENSIONS
            else {
                modelNames <- .Mclust$emModelNames
                if (n <= d) {
                  m             <- match(c("EEE", "EEV", "VEV", "VVV"), .Mclust$emModelNames,nomatch = 0)
                  modelNames    <- modelNames[-m]
                }
            }
        }
        # NUMBER OF CLUSTERS TO LOOK FOR, MAY BE A VECTOR OF NUMBERS
        if (is.null(G)) {
            # IF THE MODELING WAS CHOSEN WITHOUT POISSON NOISE, THEN DEFAULT G
            # RANGE FROM 1 TO 9, ELSE 0 IS ALLOWED FOR NOISE MODELING ONLY (SEE
            # DOCUMENTATION)
            G           <- if (is.null(initialization$noise)) 
                1:9
            else 0:9
        }
        else {
            G           <- sort(as.numeric(G))
        }
        # ORDERED CLUSTER NUMBERS AND MODELS
        Gall            <- G
        Mall            <- modelNames
    }
    #
    # CONTROL PROVIDED PARAMETERS
    #
    else {
        Glabels         <- dimnames(x)[[1]]
        Mlabels         <- dimnames(x)[[2]]
        if (is.null(G)) 
            G           <- Glabels
        if (is.null(modelNames)) 
            modelNames  <- Mlabels
        # 
        # ARE THE PROVIDED PARAMETERS POSSIBLE, MATCH CONTROL
        #
        Gmatch                  <- match(as.character(G), Glabels, nomatch = 0)
        Mmatch                  <- match(modelNames, Mlabels, nomatch = 0)
        if (all(Gmatch) && all(Mmatch)) {
            attr(x, "G")                <- as.numeric(G)
            attr(x, "modelNames")       <- modelNames
            attr(x, "returnCodes")      <- attr(x, "returnCodes")[as.character(G), modelNames, drop = FALSE]
            return(x[as.character(G), modelNames, drop = FALSE])
        }
        Gall            <- sort(as.numeric(unique(c(as.character(G), Glabels))))
        Mall            <- unique(c(modelNames, Mlabels))
    }
    if (any(as.numeric(G)) < 0) {
        if (is.null(initialization$noise)) {
            stop("G must be positive")
        }
        else {
            stop("G must be nonnegative")
        }
    }
    #
    # IF DIM OF THE DATA IS ONE, THEN A MATCHING ON THE PROVIDED MODELS TO
    # COMPUTE ON IS DONE E, V
    #
    if (d == 1 && any(nchar(modelNames) > 1)) {
        Emodel          <- any(sapply(modelNames, function(x) charmatch("E", x, nomatch = 0)[1]) == 1)
        Vmodel          <- any(sapply(modelNames, function(x) charmatch("V", x, nomatch = 0)[1]) == 1)
        modelNames      <- c("E", "V")[c(Emodel, Vmodel)]
    }
    #
    # INITIALIZE VARIABLES
    #
    l                   <- length(Gall)
    m                   <- length(Mall)
    EMPTY               <- -.Machine$double.xmax
    BIC  <- AIC  <- SIG <- RET <- matrix(EMPTY, nrow = l, ncol = m, dimnames = list(as.character(Gall), as.character(Mall)))
    if (!is.null(x)) {
        AIC[dimnames(x)[[1]], dimnames(x)[[2]]]         <- x
        BIC[dimnames(x)[[1]], dimnames(x)[[2]]]         <- x
        SIG[dimnames(x)[[1]], dimnames(x)[[2]]]         <- x
        RET[dimnames(x)[[1]], dimnames(x)[[2]]]         <- attr(x, "returnCodes")
        AIC                                             <- AIC[as.character(G), modelNames, drop = FALSE]
        BIC                                             <- BIC[as.character(G), modelNames, drop = FALSE]
        RET                                             <- RET[as.character(G), modelNames, drop = FALSE]
        SIG                                             <- SIG[as.character(G), modelNames, drop = FALSE]
    }
    G                   <- as.numeric(G)
    Glabels             <- as.character(G)
    Gout                <- G
    stats				<- NULL
    #
    # LIST TO STORE THE RESULTS OF THE MULTIVARIATE NORMAL FITS
    #
    outputList          <- NULL
    #
    if (is.null(initialization$noise)) {
    	#
    	# HERE!!!
    	#
        if (G[1] == 1) {
            for (mdl in modelNames[BIC["1", ] == EMPTY]) {
                out             <- mvn(modelName = mdl, data = data, prior = prior)
                out.scores      <- mclust2_bicaic(modelName = mdl, loglik = out$loglik, n = n, d = d, G = 1, equalPro = FALSE)
                BIC["1", mdl]   <- out.scores[1] 
                AIC["1", mdl]   <- out.scores[2] 
                RET["1", mdl]   <- attr(out, "returnCode")
                outputList    	<- c(outputList,list(out))                
            }
            if (l == 1) {
                BIC[BIC == EMPTY] <- NA
#                return(structure(BIC, BIC=BIC, AIC=AIC,  SIG=SIG, NSignificant = NSignificant, G = G, 
                return(structure(BIC, BIC=BIC, AIC=AIC,  SIG=SIG, G = G, 
                				modelNames = modelNames, prior = prior, control = control, 
                                initialization = list(hcPairs = initialization$hcPairs, subset = initialization$subset), 
                                warn = warn, n = n, d = d, oneD = oneD, returnCodes = RET, out = outputList, 
                                class = "mclust2_BIC"))
            }
            G                   <- G[-1]
            Glabels             <- Glabels[-1]
        }
        if (is.null(initialization$subset)) {
            if (is.null(initialization$hcPairs)) {
                if (d != 1) {
                  if (n > d) {
                    hcPairs     <- hc(modelName = "VVV", data = data)
                  }
                  else {
                    hcPairs     <- hc(modelName = "EII", data = data)
                  }
                }
                else {
                  hcPairs       <- NULL
                }
            }
            else hcPairs        <- initialization$hcPairs
            if (d > 1 || !is.null(hcPairs)) 
                clss            <- hclass(hcPairs, G)
            #
            # AFTER FITTING A SINGLE MULTIVARIATE MODEL ON THE DATA (1-CLUSTER, PREVIOUS LOOP),
            # CALCULATE BY EM ALGORITHM, EVENTUALY INITIALIZED BY HIERARCHICAL CLUSTERING, THE 
            # STATISTIC FOR "G" CLUSTERS IN THE DATA.
            #
            for (g in Glabels) {
                if (d > 1 || !is.null(hcPairs)) {
                  z             <- unmap(clss[, g])
                }
                else {
                  z             <- unmap(qclass(data, as.numeric(g)))
                }
                for (modelName in modelNames[BIC[g, ] == EMPTY]) {
                  out           <- me(modelName = modelName, data = data, z = z, prior = prior, control = control, warn = warn)
                  out.scores    <- mclust2_bicaic(modelName = modelName, loglik = out$loglik, n = n, d = d, G = as.numeric(g), equalPro = control$equalPro)
                  BIC[g, modelName] <- out.scores[1] 
                  AIC[g, modelName] <- out.scores[2] 
                  RET[g, modelName] <- attr(out, "returnCode")
                  if(RET[g, modelName] == 0){
                  	stats <- NULL
               		# stats 		<- mclust2_statistics(mclust2_mergeSibpairs(matrix(map(out$z),
                  	#							dimnames=list(dimnames(out$z)[[1]], NULL))))
             	    # SIG[g, modelName]	<- stats["Significant?","T"]
                  }
                  else{
                  	stats				<- NULL
                  	SIG[g, modelName]	<- 0
                  	}
                    
                  outputList    <- c(outputList,list(out))
                }
            }
        }
        else {
            if (is.logical(initialization$subset)) 
                initialization$subset <- (1:n)[initialization$subset]
            if (is.null(initialization$hcPairs)) {
                if (d != 1) {
                  if (n > d) {
                    hcPairs     <- hc(modelName = "VVV", data = data[initialization$subset, ])
                  }
                  else {
                    hcPairs     <- hc(modelName = "EII", data = data[initialization$subset, ])
                  }
                }
                else {
                  hcPairs <- NULL
                }
            }
            else hcPairs <- initialization$hcPairs
            if (d > 1 || !is.null(hcPairs)) 
                clss <- hclass(hcPairs, G)
            for (g in Glabels) {
                if (d > 1 || !is.null(hcPairs)) {
                  z <- unmap(clss[, g])
                }
                else {
                  z <- unmap(qclass(data[initialization$subset], as.numeric(g)))
                }
                dimnames(z) <- list(as.character(initialization$subset), NULL)
                for (modelName in modelNames[!is.na(BIC[g, ])]) {
                  ms 			<- mstep(modelName = modelName, z = z, data = as.matrix(data)[initialization$subset, ], 
                  						prior = prior, control = control, warn = warn)
                  es 			<- do.call("estep", c(list(data = data, warn = warn), ms))
                  out 			<- me(modelName = modelName, data = data, z = es$z, prior = prior, control = control, warn = warn)
                  out.scores    <- mclust2_bicaic(modelName = modelName, loglik = out$loglik, n = n, d = d, G = as.numeric(g), 
                  						equalPro = control$equalPro)
                  BIC[g, modelName] <- out.scores[1] 
                  AIC[g, modelName] <- out.scores[2] 
                  RET[g, modelName] <- attr(out, "returnCode")
                  if(RET[g, modelName] == 0){
                  		stats 	<- NULL
#               		 stats 		<- mclust2_statistics(mclust2_mergeSibpairs(data,matrix(map(out$z),
#                  								dimnames=list(dimnames(out$z)[[1]], NULL))))
             	     # SIG[g, modelName]	<- stats["Significant?","T"]
                  }
                  else{
                  	stats			<- NULL
                  	SIG[g, modelName]	<- 0
                  	}
                  outputList    <- c(outputList,list(out))
                }
            }
        }
    }
    structure(BIC, BIC=BIC, AIC=AIC, SIG=SIG, G = Gout, modelNames = modelNames, prior = prior, control = control, 
        initialization = list(hcPairs = hcPairs, subset = initialization$subset, noise = initialization$noise), 
        Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD, 
        returnCodes = RET, out = outputList,class = "mclust2_BIC")
}

`mclust2_bicaic` <-
function (modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE,...)
{
   modelName    <- switch(EXPR = modelName, XII = "EII", XXI = "EEI",XXX = "EEE", modelName)
   if (G == 0) {
      if (!noise)
      stop("undefined model")
      nparams <- 1
   }
   else {
         nparams <- nVarParams(modelName, d, G) + G * d
      if (!equalPro)
         nparams <- nparams + (G - 1)
      if (noise)
         nparams <- nparams + 2
   }
   #
   # RETURN VECTOR WITH LOG BIC AND AIC SCORES
   #
   return(c(2 * loglik - nparams * logb(n), 2 * loglik - 2 * nparams))
}

`model_based_clustering` <-
function(data, config  # em_clustering OR mclust2
   , G_set , model_names=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"), clust_fun = mclust2) {
   prefix <- attr(config,"prefix")
   # DATA FOR THE CLUSTERING (CC AND TRANSFORMED)
   canalysis_data        <- get_canalysis_data(data)
   # CLUSTERING
   cresult	        <- clust_fun(data=canalysis_data,G=G_set,modelNames=model_names)
   # BIC TABLE AND EVENTUALLy OUTPUT INTO A CSV
   cresult[["bic_table"]]	        <- attr(cresult[["out"]],"BIC")
   cresult[["bic_table_relative"]]     <- 100*(cresult[["bic_table"]]/max(cresult[["bic_table"]],na.rm=TRUE)-1)
   write.table(file=paste2(prefix,"_Relative_BIC_Table.csv"),cresult[["bic_table_relative"]],dec=",",sep=";")
   # SIMPLIFY MCLUST2 DATA STRUCTURE
   # the two fun_clust do not produce exactly the same data structure ......
   tmp_cresult_out <- attr(cresult[["out"]],"out")
   if(!is.null(tmp_cresult_out))
      cresult[["out"]] <- tmp_cresult_out
   return(cresult)
}

`model_knn` <-
function(formula, trainset, testset, k = 1, l = 0, prob = FALSE, use.all = TRUE){
   return(knn(trainset[,-match("class",names(trainset))], testset[,-match("class",names(testset))] ,trainset[,"class"], k = k, l = l, prob = prob, use.all = use.all))
}

`model_naive_bayes` <-
function(formula, trainset, testset, laplace = 1){
   return(naiveBayes(class ~ . , data=trainset, laplace = 1))
}

`model_svm` <-
function(formula, trainset, testset, kernel="linear",type="C-classification"){
   return(svm(class ~ . , data=trainset,kernel=kernel,type=type))
}

`paste2` <-
function(...){
	return(paste(...,collapse="",sep=""))
}

`plot_init` <-
function(visu_params,title_str){
   plot(x=0 ,new=TRUE ,ann=FALSE ,pch=18,col="white",axes=FALSE
      , xlim=c(visu_params[["xlim"]]$min,visu_params[["xlim"]]$max)
      , ylim=c(visu_params[["ylim"]]$min,visu_params[["ylim"]]$max))
   title(main=title_str)
}

`plot_legend` <-
function(pattern, color_settings, xcoord, ycoord){
   pattern <- pattern[[1]]
   for(g in 1:nrow(pattern)){
      g_name <- color_settings[["cluster_order"]][g]
      y1 <- y0 <-  ycoord +  5 * g / nrow(pattern)
      x0 <-        xcoord
      x1 <-        1.2 * xcoord
      arrows(x0,y0,x1,y1,col=color_settings[["cluster_color"]][g],length=0,lwd=3)
      text(x0+0.15,y1,labels=paste2(g_name," (",pattern[as.numeric(g_name)],")"),pos=4)
   }
}

`plot_parcoord` <-
function(pattern,T_s,color_sel, lty="full", lwd=3,title_prefix = NULL){
   # LOOP OVER EACH GROUP PATTERN THAT IS TO VISUALIZE
   for(g in color_sel[["cluster_order"]]){
      # LOOP OVER THE DIFFERENT STATISTICS FOR EACH PATTERN
      D_i_s             <- cbind(X=as.numeric(pattern[g,row.names(T_s)]),Y=as.numeric(T_s[,"visu_ycoord"]))
      D_i_s	        <- D_i_s[sort.list(D_i_s[,"Y"]),]
      gap               <- D_i_s[2,"Y"] - D_i_s[1,"Y"]
      for(l in 1:(nrow(D_i_s)-1))
         if(!(D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap))
            arrows(D_i_s[l,1], D_i_s[l,2], D_i_s[l+1,1], D_i_s[l+1,2],
               col=color_sel[["cluster_color"]][g], length=0, lwd=lwd,
               lty=lty)
         # ELSE, AS (is_white_gap == TRUE) THEN DO NOT DRAW ANY ARROW...
   }
}

`predict_knn` <-
function(model, testset){
   return(model)
}

`predict_naive_bayes` <-
function(model,testset){
   return(map(predict(model,testset[,-match("class",names(testset))],type="raw")))
}

`predict_svm` <-
function(model, testset){
   predict(model,testset[,-match("class",names(testset))])
}

`retrieve_model` <-
function(canalysis,query){
   # init
   params       <- list()
   # retrieve the model corresponding to the query
   mbc_params   <- strsplit(query,",")[[1]]
   params       <- list(model=mbc_params[1],G=as.numeric(mbc_params[2]))
   i            <- 1 
   while(!(canalysis[["out"]][[i]]$G == params$G & canalysis[["out"]][[i]]$modelName == params$model & i <= length(canalysis[["out"]])))
      i         <- i+1
   params[["model"]]     <- canalysis[["out"]][[i]]
   params[["labelling"]] <- map(params[["model"]]$z,warn=FALSE)
   return(params)
}

`save_data_rep` <-
function(version = 1.3){
   fileName <- paste2("PKG_",version,".RData")
   save(file = fileName, list = c("data_OA", "data_OA_bootstrap_idx",
      "data_PD_4y", "data_PD_y_1", "data_PD_y_2", "data_PD_y_3", "data_PD_y_4",
      "data_PHARMAIT", "setup_OA_data", "setup_OA_visu_params", "setup_PD_data",
      "setup_PD_visu_params"))
}

`save_model` <-
function(model, prefix){
   data_out        <- cbind(  model[["z"]]
                            , class=map(model[["z"]]))
   colnames(data_out)[1:(ncol(data_out)-1)] <- 1:(ncol(data_out)-1)
   fileNameCSV     <- paste2(prefix,"_Class_",model[["modelName"]],"-",model[["G"]],".csv") 
   print(fileNameCSV)
   write.csv2(data_out,file=fileNameCSV)
}

`set_calgo` <-
function(calgo = "model_based_clustering"
        , args = list( G_set = 1:9 , model_names=c("EII", "VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"))){
   clust_fun <- function(data=NULL,config=NULL){
      return(do.call(calgo,append(list(data=data,config=config),args)))
      }
   return(clust_fun)
}

`set_canalysis_config` <-
function(
          prefix = NULL
        , visu_params = NULL
        , stats_fun = list(oddratios = get_stats_fun(fun_name="oddratios"))
        , clustering_fun = model_based_clustering
        , bbox_threshold = 5){
   #
   prefix <- paste2(today(),"_",prefix)
   #
   cdata_out <- structure(
     visu_params
   ,  prefix = prefix
   , stats_fun = stats_fun
   , clustering_fun = clustering_fun
   , bbox_threshold = bbox_threshold
   , class = "canalysis_config"
   )
   return(cdata_out)
}

`set_canalysis_data` <-
function(data, data_o = NULL, tdata = NULL, settings, init_fun = list(init_data_cc)){
   tdata_list   <- list()
   settings     <- as.matrix(settings)
   # BACKUP DATA INTO DATA_O 
   if(is.null(data_o) )
      data_o <- data
   # DO INIT FUNCTIONS OF THE DATA
   for(fun in init_fun)
      data <- fun(data,settings)
   # EVENTUALLY TRANSFORM THE DATA
   if(is.null(tdata)){
      t_result  <- transform_canalysis_data(data,settings)
      data      <- t_result[["data"]]
      tdata     <- t_result[["tdata"]]
   }
   #
   cdata_out <- structure( data , data_o = data_o  , tdata = tdata , 
                        settings = settings , class = "canalysis_data")
   return(cdata_out)
}

`set_canalysis_result` <-
function(data = NULL, config = NULL , result = list()){
   if(class(data) == "canalysis_data" & class(config) == "canalysis_config"){
      cdata_out <- structure(
      result
      , config = 
      , data = data
      , model = model
      , settings = settings
      , class = "canalysis_result"
      )
   }
   else
      cdata_out <- NULL
   return(cdata_out)
}

`set_plot_data` <-
function(
           canalysis_data 
         , canalysis_config 
         , model = NULL
         , plot_data_fun = list(center_pattern  = mean)
#                            , lowquant         = function(x){return(quantile(x,probs = seq(0, 1, 0.025))[1])}
#                            , upquant          = function(x){return(quantile(x,probs = seq(0, 1, 0.025))[40])}
#                            , center_pattern   = median
         , title_str = NULL
         , type = "plot_parcoord"){
   # RETRIEVE GLOBAL SETTINGS, AND VISUALIZATION PARAMETERS
   settings             <- attr(canalysis_data,"settings")
   # MAKE LABELLED DATA_SET
   labelled_data        <- as.matrix(cbind(canalysis_data,class=map(model[["z"]])))
   # DETERMINES ROW AND COLUMN ORDERINGS OF THE 'CENTER' PATTERN
   center_pattern       <- compute_plot_data(labelled_data,plot_data_fun[[1]])
   cluster_dendro       <- hclust(dist(center_pattern))
   var_dendro           <- hclust(dist(t(center_pattern)))
   center_pattern       <- center_pattern[cluster_dendro$order, rev(sort(colnames(center_pattern)))]
   # SET PARCOORD AND LEGEND COLOR SETTINGS 
   color_settings       <- get_coloring_scheme(labelled_data)
   # TITLE:     DEFINE VISUALIZATION TITLE 
   # 
   if(type == "plot_parcoord"){
      # FOR PARCOORD, LIMIT THE VISU TO THE VAR OF SUCH A VIGNETTE
      visu_var          <- row.names(settings[settings[,"var_by_groups"] == title_str,])
      # COMPUTE THE DATA_PLOT
      plot_data            <- list()
      title_str         <- paste2("(",model[["modelName"]],",",model[["G"]],") ", title_str)
      for(fun_name in names(plot_data_fun))
         plot_data[[fun_name]] <- compute_plot_data(as.matrix(labelled_data[,c(visu_var,"class")]),plot_data_fun[[fun_name]])
      plot_fun <- function(x, visu_params = canalysis_config, color_sel = color_settings, T_s = settings[visu_var,],title_prefix = title_str){
                        plot_init(visu_params,title_prefix)
                        # LOOP OVER THE DIFFERENT STATS TO PLOT == RAYS TO GRAPH
                        for(x_name in names(x)){
                           pattern         <- x[[x_name]][color_sel[["cluster_order"]],]
                           lwd             <- 3
                           lty             <- "solid"
                           if(x_name != "center_pattern") {
                              lwd          <- 1
                              lty          <- "dashed"
                           }
                           plot_parcoord(pattern = pattern, T_s = T_s,
                              color_sel = color_sel, lwd = lwd, lty = lty,
                              title_prefix = title_prefix)
                        }
                        axis(2,at=as.numeric(T_s[,"visu_ycoord"]),labels=colnames(pattern),las=2,tick=FALSE)
                        axis(1,at=seq(from = visu_params[["xlim"]]$min , to =
                           visu_params[["xlim"]]$max , by = ((visu_params[["xlim"]]$max -
                           visu_params[["xlim"]]$min)/4)))
                  }
   }
   if(type == "plot_image"){
      plot_data <- list(center_pattern=center_pattern)
      plot_fun  <- function(x, visu_params = canalysis_config,title_prefix = title_str){
                        plot_init(visu_params,title_prefix)
                        x <- x[[1]]
                        image(1L:nrow(x), 1L:ncol(x), x, xlim = 0.5 + c(0, nrow(x)),
                           ylim = 0.5 + c(0, ncol(x)), axes = FALSE, xlab = "",
                           ylab = "", col=visu_params[["color_gradient"]],
                           main="Average pattern visualization")
                        axis(2, 1L:ncol(x), labels = colnames(x), las = 2,
                           line = -0.5, tick = 0, cex.axis = visu_params[["cex"]])
                        axis(1, 1L:nrow(x), labels = row.names(x), las =
                           1, line = -0.5, tick = 0, cex.axis = visu_params[["cex"]])
                        }
   }
   if(type == "plot_legend"){
      plot_data <- list(class_count = table(labelled_data[,"class"]))
      plot_fun  <- function(x, visu_params = canalysis_config){
                        plot_legend(pattern = x, color_settings =
                        color_settings, xcoord= canalysis_config[["legend_xcoord"]], 
                        ycoord= canalysis_config[["legend_ycoord"]])
                     }
   }
   if(type == "plot_dendro_cluster"){
      plot_data <- list(cluster_dendro = cluster_dendro)
      plot_fun  <- function(x, visu_params = canalysis_config,title_prefix = title_str){
                        plot_init(visu_params,title_prefix)
                        plot(x[[1]], axes = FALSE, yaxs = "i", main="Clusters arranged
                        by similarity", ylab=NULL, cex=visu_params[["cex"]])
                     }
   }
   if(type == "plot_dendro_var"){
      plot_data <- list(var_dendro = var_dendro)
      plot_fun  <- function(x, visu_params = canalysis_config,title_prefix = title_str){
                        plot_init(visu_params,title_prefix)
                        plot(x[[1]], axes = FALSE, yaxs = "i", main="Clusters arranged
                        by similarity", ylab=NULL, cex=visu_params[["cex"]])
                     }
   }
   # CONSTRUCT THE DATA STRUCTURE
   plot_data_out <- structure(plot_data, plot_fun = plot_fun, class = "plot_data")
   return(plot_data_out)
}

`set_tdata` <-
function(testimate=NULL,tfun=NULL){
  tdata_out <- structure(
       testimate
     , tfun  = tfun
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

`stability_assessment` <-
function(
     data               = na.omit(as.data.frame(df))
   , canalysis_variables  = NULL
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
                                                transform_addnoise(data, model, canalysis_variables = canalysis_variables, type, 
                                                relative_noise_amount = noise_vector[noise_idx], rand_seed = 6013+3*i))}
         filePrefix_local <- paste2(filePrefix,noise_vector[noise_idx],"_",adj,"_")
         # PROCEED TO CLUSTER ANALYSIS 10 TIMES WITH A DIFFERENT SEED
         for(i in 1:10)
            canalysis_results[[adj]][[noise_idx]][[i]] <- analysis(
                                                              data                   = data
                                                            , canalysis_variables      = canalysis_variables 
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
   quantile_df <- matrix(NA,length(noise_vector),11,dimnames=list(100*noise_vector,c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%")))
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
      axis(1,at=seq(0,1,0.11),labels=sprintf("%.1f%%",as.numeric(row.names(quantile_a))),las=2,cex.axis=1.7)
      axis(2,at=seq(0,1,0.1),labels=colnames(quantile_a),cex.axis=1.7)
      }
   graphics.off()
   return(list(canalysis=canalysis_results,quantile_a))
}

`statistical_patterns` <-
function(data,fun_list=list(avg=mean,median=median)){
   data_out <- list()
   for(fun in names(fun_list))
      data_out[[fun]] <- compute_plot_data(data,fun=fun_list[[fun]])
   return(data_out)
}

`statistics_lambdasibs` <-
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

`statistics_logodds` <-
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

`stratified_cv` <-
function(data,K=10){
   indexes <- list()
   G <- length(table(data$class))
   g <- as.numeric(names(sort(table(data$class))[1]))
   row_number <- nrow(data[data$class == g,])
   out_cv <- stratified_cv_folds(row_number,K=K,g)
   while(is.null(out_cv)){
      K <- K-1
      print(K)
      out_cv <- stratified_cv_folds(row_number,K=K,g)
   }
   for(g in 1:G)
      indexes[[g]] <- stratified_cv_folds(row_number,K,g)
   return(list(indexes=indexes,K=K))
}

`stratified_cv_folds` <-
function(n, K, G) 
{
   set.seed(6013)
   # WE DO A ROTATION ON THE ROWS OF THE FINAL FOLD MATRIX FOR NOT HAVING
   # ALWAYS THE SAME FOLD-STRATAS WITH ONE-LESS COUNT, I.E. THE ONE FOLDS AT
   # THE END OF THE MATRIX, 
   if(G == K) 
      G <- 1
   rowPermutation <- c(((G %% K)+1):K,1:(G %% K))
   # 
   size <- n/K
   cv <- matrix(0, K, ceiling(size))
   if (size < 5 & size != 1){
      print("The number of folds is too large for the data size")
      return(NULL)
   }
   if (size == 1) 
      return(matrix(1:n, nrow = n))
   size.int <- floor(size)
   size.vector <- rep(size.int, K)
   if (size.int != size) 
      size.vector[1:((size - size.int) * K)] <- size.vector[1:((size - size.int) * K)] + 1
   group.index <- c()
   for (j in 1:K) 
      group.index <- c(group.index, rep(j, size.vector[j]))
   group.index <- group.index[sample(n, n, replace = FALSE)]
   for (j in 1:K) {
      whichj <- which(group.index == j)
      cv[j, 1:length(whichj)] <- whichj
   }
   return(cv[rowPermutation,])
}

`stratified_traintest_split` <-
function(data,K=10,perc=0.7){
   indexes <- list()
   G <- length(table(data$class))
   for(g in 1:G){
      strata_idx <- nrow(data[data$class == g,])
      indexes[[g]] <- matrix(0,0,strata_idx * (1-perc)) 
      set.seed(6013)
      for(k in 1:K)
            indexes[[g]] <- rbind(indexes[[g]],sample(1:strata_idx)[1:(strata_idx * (1-perc))])
      }
   return(list(indexes=indexes,K=K))
}

`summary_mclust2_BIC` <-
function (object, data, G = NULL, modelNames = NULL, ...) 
{
    dimData             <- dim(data)
    oneD                <- is.null(dimData) || length(dimData[dimData > 1]) == 
        1
    if (!oneD && length(dimData) != 2) 
        stop("data must be a vector or a matrix")
    if (oneD) {
        data            <- drop(as.matrix(data))
        n               <- length(data)
        d               <- 1
    }
    else {
        data            <- as.matrix(data)
        n               <- nrow(data)
        d               <- ncol(data)
    }
    initialization      <- attr(object, "initialization")
    hcPairs             <- initialization$hcPairs
    subset              <- initialization$subset
    prior               <- attr(object, "prior")
    control             <- attr(object, "control")
    warn                <- attr(object, "warn")
    oldClass(object)    <- NULL
    attr(object, "prior")               <- attr(object, "warn") <- NULL
    attr(object, "modelNames")          <- attr(object, "oneD") <- NULL
    attr(object, "initialization")      <- attr(object, "control") <- NULL
    d                   <-      if (is.null(dim(data))) 1
                                else ncol(data)
    if (is.null(G)) 
        G               <- dimnames(object)[[1]]
    if (is.null(modelNames)) 
        modelNames      <- dimnames(object)[[2]]
    bestBICs            <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
    if(is.na(bestBICs)[1])
    	return(NULL)
    temp                <- unlist(strsplit(names(bestBICs)[1], ","))
    bestModel           <- temp[1]
    G                   <- as.numeric(temp[2])
    if (G == 1) {
        out             <- mvn(modelName = bestModel, data = data, prior = prior)
        ans             <- c(list(bic = bestBICs, classification = rep(1, n), uncertainty = rep(0, n)), out)
        orderedNames    <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty")
        return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                prior = prior, control = control, initialization = initialization, class = "summary_mclust2_BIC"))
    }
    if (is.null(subset)) {
        if (d > 1 || !is.null(hcPairs)) {
            z           <- unmap(hclass(hcPairs, G))
        }
        else {
            z           <- unmap(qclass(data, G))
        }
        out             <- me(modelName = bestModel, data = data, z = z, prior = prior, control = control, warn = warn)
    }
    else {
        if (d > 1 || !is.null(hcPairs)) {
            z           <- unmap(hclass(hcPairs, G))
        }
        else {
            z           <- unmap(qclass(data[subset], G))
        }
        ms              <- mstep(modelName = bestModel, prior = prior, z = z, data = as.matrix(data)[subset, ], control = control, warn = warn)
        es              <- do.call("estep", c(list(data = data), ms))
        out             <- me(modelName = bestModel, data = data, z = es$z, prior = prior, control = control, warn = warn)
    }
    obsNames <- if (is.null(dim(data))) {
        names(data)
    }
    else {
        dimnames(data)[[1]]
    }
    classification      <- map(out$z)
    uncertainty         <- 1 - apply(out$z, 1, max)
    names(classification) <- names(uncertainty) <- obsNames
    ans                 <- c(list(bic = as.vector(bestBICs[1]), classification = classification, uncertainty = uncertainty), out)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", "z", "classification", "uncertainty")
    structure(ans[orderedNames], bestBICvalues = bestBICs, prior = prior, 
        control = control, initialization = initialization, class = "summary_mclust2_BIC")
}

`today` <-
function(){
   tmp_date <- format(Sys.time(), "%Y-%m-%d")
   return(tmp_date)
}

`transform_ABSMAX` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun         <- list(function(x){return(max(abs(x,na.rm=TRUE),na.rm=TRUE))}
                      , function(x){return(do.call('/',x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_ALL` <-
function(vdata,tfun=list(function(x){return(0)},function(x){return(do.call('+',x))}),tdata){
   if(length(tdata) == 0) 
      testimate         <- tfun[[1]](vdata)
   else 
      testimate         <- as.numeric(tdata)
   vdata                <- tfun[[2]](list(vdata,testimate))
   tdata_out            <- set_tdata(testimate=testimate,tfun=tfun)
   return(list(data=vdata,tdata=tdata_out))
}

`transform_AVG` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun         <- list(function(x){return(mean(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_L1` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun      <- list(function(x){return(x/sqrt(sum(x,na.rm=TRUE)))}
                      , function(x){return(do.call('/',x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_L2` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun      <- list(function(x){return(x/sqrt(sum(x^2,na.rm=TRUE)))}
                      , function(x){return(do.call('/',x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_MAX` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun         <- list(function(x){return(max(x,na.rm=TRUE))}
                      , function(x){return(do.call('-', x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_MEDIAN` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun         <- list(function(x){return(median(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_MIN` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun       <- list(function(x){return(min(x,na.rm=TRUE))}
                      , function(x){return(do.call('-',x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_SIGMA` <-
function(data,tdata,var){
   vdata        <- data[,var]
   tfun         <- list(function(x){return(sd(x,na.rm=TRUE))}
                      , function(x){return(do.call('/', x))})
   return(transform_ALL(vdata=vdata,tfun=tfun,tdata=tdata))
}

`transform_addnoise` <-
function(
        data,
        model=NULL,
        canalysis_variables = NULL,
        tfun="addnoise",
        relative_noise_amount = 0.01,
        rand_seed = 6013
        )
{
   if(is.null(model)){
      for(o in which(canalysis_variables %in% colnames(data))){
         local_name             <- colnames(data)[o]
         local_rand_seed        <- rand_seed+o
         set.seed(local_rand_seed)
         local_rand_vector      <- data.frame(rnorm(    length(data[,o])
                                                        , mean = mean(data[,o])
                                                        , sd = relative_noise_amount * sd(data[,o])
                                                   )
                                             , row.names = row.names(data))
         colnames(local_rand_vector) <- "values"
         model[[local_name]]    <- list(  rand_seed = local_rand_seed
                                        , rand_vector = local_rand_vector
                                        , relative_noise_amount = relative_noise_amount)
         # (EXPECT NO NA'S WHILE SUMMING BECAUSE "0.54345<-NA+0.54345"...)
         data[,o] <- data[,o] + local_rand_vector[,1]
      }
   }
   else{
      for(o in colnames(data))
         data[,o] <- data[,o] + model[[o]][["rand_vector"]][row.names(data),]
   }
   return(list(data=data,tfun=tfun,model=model))
}

`transform_adjust` <-
function(data, tdata, tformula){
   var          <- strsplit(tformula,"~")[[1]][1]
   data         <- as.data.frame(data)
   if(is.null(tdata)){
      tdata             <- lm(as.formula(tformula),data=data)
      tdata[["formula"]]<- tformula
      tdata[["formula"]]<- tformula
      vdata             <- residuals(tdata)
   }
   else
      vdata             <- data[,var]-predict(tdata,data)
   return(list(data=vdata,tdata=tdata))
}

`transform_canalysis_data` <-
function(data,settings){
   tdata_out <- list()
   # FOR EACH VAR
   for(var in row.names(settings)){
      fun_transform     <- strsplit(settings[var,"fun_transform"],"[ ]+")[[1]]
      tdata_out[[var]]  <- list()
      # TODO transform the char string into a function "Call", and test
      for(lfun in fun_transform){
         pattern        <-paste2("[('\")]")
         lfun           <- strsplit(gsub(pattern," ",lfun),"[ ]+")[[1]]
         arg_list       <- list(data,tdata=attr(data,"tdata")[[var]])
         if(length(lfun) == 1) 
            lfun        <-   append(lfun,var)
         arg_list       <- append(arg_list,lfun[2])
         # RETURNS 'VDATA' AND TDATA_OUT
         var_result     <- do.call(lfun[1],args=arg_list)
         # UPDATE THE MAIN DATA FRAME
         data[,var]     <- as.numeric(var_result[["data"]][row.names(data),])
         # KEEP RECORD OF EACH INDIVIDUAL TRANSFORMATION FOR THE 'VAR'
         tdata_out[[var]][[lfun[1]]]     <- var_result[["tdata"]]
      }
   }
   return(list(data=data,tdata=tdata_out))
}

`transform_longitudinal` <-
function(data, fun, canalysis_variables = parkinson_canalysis_variables,effect="duration"){
   data_out             <- data
   m                    <- fun(as.data.frame(data[,c(canalysis_variables,effect),1]))
   m_effect             <- 
   data_out[row.names(m[["data"]]),colnames(m[["data"]]),1]     <- data.matrix(m[["data"]])
   for(i in 2:dim(data)[[3]]){
      # WE WANT TO REMOVE THE EFFECT OF DISEASE DURATION (AKA THE TIME) FOR BASELINE,
      # BUT WE PRESERVE THE TIME EFFECT IN THE LONGITUDINAL ANALYSIS. 
      data_tmp          <- as.data.frame(cbind(data[,canalysis_variables,i],data[,effect,1]))
      colnames(data_tmp)[ncol(data_tmp)] <- effect 
      m_tmp          <- fun(data_tmp,model = m[["model"]])
      data_out[row.names(m_tmp[["data"]]),canalysis_variables,i]  <- data.matrix(m_tmp[["data"]])[,canalysis_variables]
      }
   return(data_out)
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

`variable_graphic_histograms` <-
function(data, which_data = NULL ,prefix=NULL){
   # GET CANALYSIS DATA
   ldata        <- get_canalysis_data(data, which_data = which_data)
   ldata        <- as.data.frame(ldata)
   hist_list    <- hist_stats <- list()
   hist_counts  <- c()
   models       <- attr(data,"model")
   # INIT PLOTTING REGION
   par('mai'=c(0.5,0.5,0.5,1),'mfrow'=c(4,5))
   # LOOP OVER VARIABLES SELECTED FOR THE CANALYSIS
   for(i in colnames(ldata))
   {
      hist_list[[i]]  <- hist(ldata[,i],plot=FALSE)
      local_stat_txt <- ""
      hist_stats[[i]] <- list(txt="",xlim=list(),ylim=list())
      for(j in names(models)){
         local_model <- models[[j]][["model"]][[i]]
         if(models[[j]][["type"]] == "lm")
            local_stat_txt <- paste2(local_stat_txt,local_model[["formula"]],"\n")
         else
            local_stat_txt <- paste2(local_stat_txt,sprintf("%1.2f",local_model[[1]]),", ")
      }
      local_stat_txt <- paste2(local_stat_txt,length(table(ldata[,i]))," distinct values")
      hist_stats[[i]][["txt"]] <- local_stat_txt
      hist_stats[[i]][["xlim"]] <- c(min(ldata[,i],na.rm=TRUE),max(ldata[,i],na.rm=TRUE))
#           hist_stats[[i]][["xlim"]] <- c(quantile(data.matrix(ldata[,i]),probs=c(0.05),na.rm=TRUE),
#                quantile(data.matrix(ldata[,i]),probs=c(0.95),na.rm=TRUE))
      hist_stats[[i]][["ylim"]] <- range(hist_list[[i]][["counts"]])
   }
   #
   for(i in 1:ncol(ldata))
   {
      par('lab'=c(3,3,3),  # number of ticks
            'mar'=c(2.5,2.5,1.5,1) # margins bottom,left,top,right
            )
      plot(hist_list[[i]],
            main=dimnames(ldata)[[2]][i],
            xlab="",ylab="",xlim=hist_stats[[i]][["xlim"]],ylim=hist_stats[[i]][["ylim"]]
            #,ylog=TRUE
            ) #!!!!
      text(x=mean(hist_stats[[i]][["xlim"]]),y=hist_stats[[i]][["ylim"]][2]*0.85,hist_stats[[i]][["txt"]])
      }
}

`variable_graphic_report` <-
function(data,which_data=NULL,prefix=NULL){
   if(class(data) != "canalysis_data") {
      print("Error (variables_graphic_report): You must provide data of type canalysis_data")
      break ;
   }
   else {
      # CORRELATION MATRIX AND ABSOLUTE VALUED CORRELATION MATRIX TAKEN AS A
      # DISTANCE MATRIX TO PLOT A DENDROGRAM WHICH ILLUSTRATES OUTCOME'S RELATION
      ldata     <- get_canalysis_data(data, which_data = which_data)
      ldata     <- as.data.frame(ldata[1:nrow(ldata),])
      cor_mat   <- cor(ldata,use="pairwise.complete.obs")
      write.table(cor_mat,file=paste2(prefix,"Correlation_Matrix.csv"),dec=",",sep=";")
      hclust(as.dist(abs(cor_mat)))
      # BOXPLOTS SUMMARY DITRIBUTIONAL SHAPE OF THE OUTCOMES BY THEIR BOX PLOTS
      # 1ST AND 3RD QUARTILES, MEDIAN, 95% BOUNDS, + EXTREMES BEYOND THE 95%
      par(mfrow=c(1,1),'mai'=c(1.2,0.7,0.7,0.7),new=FALSE)
      boxplot(ldata,'las'=2,main=paste2(prefix,"Boxplot"))
      # PLOT HISTOGRAMS FOR EACH VARIABLE
      variable_graphic_histograms(data,which_data,prefix=prefix)
   }
}

