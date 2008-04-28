###################### BEGIN, LIBRARY & DATA IMPORTS #######################
# IMPORT LIBRARIES OF FUNCTIONS
# 3RD PARTY LIBRARIES
library(abind)
library(RColorBrewer)
library(mclust)
library(e1071)
# Return Today date
today <- function(){
   tmp_date <- format(Sys.time(), "%Y-%m-%d")
   return(tmp_date)
}
###################### END, LIBRARY & DATA IMPORTS #########################
#
###################### BEGIN, FUNCTION DEFINITION SECTION ##################
# 2008-01-13:
# 
#
save_models <- function(data,filePrefix=NULL,threshold=5){
   m_set                <- get_model_set(data[["cluster_analysis"]][["bic_table_relative"]],bbox_threshold=threshold)
   for(l in 1:length(data[["cluster_analysis"]][["out"]])){
      local_canalysis   <- data[["cluster_analysis"]][["out"]][[l]]
      g_match           <- match(local_canalysis[["G"]],m_set[,"G"])
      m_match           <- match(local_canalysis[["modelName"]],m_set[,"modelName"])
      if( !is.na(g_match) & !is.na(m_match)){
	data_out        <- cbind(         ids=row.names(data[["transformed"]])
                                        , class=map(local_canalysis[["z"]])
                                        , local_canalysis[["z"]])
        colnames(data_out)[3:ncol(data_out)] <- 1:local_canalysis[["G"]]
        fileNameCSV     <- paste2(filePrefix,"_Class_",local_canalysis[["modelName"]],"-",local_canalysis[["G"]],".csv") 
        print(fileNameCSV)
        write.csv2(data_out,file=fileNameCSV)
      }
   }
}
# 2007-12-05:
#
#
paste2 <- function(...){
	return(paste(...,collapse="",sep=""))
}
# 2007-12-10:
#
# 
transform_outcomes <- function(data,fun_transform = NULL){
   #
   # APPLY USER-PROVIDED DATA TRANSFORMATION
   #
   # data_transformed <- data[["cc"]][,unlist(data[["sumscore_groups"]])]
   data_transformed <- data[["cc"]]
   model_out <- list()
   for(fun_name in names(fun_transform)){
      model_out[[fun_name]] <- fun_transform[[fun_name]](data_transformed)
      data_transformed      <- model_out[[fun_name]][["data"]]
   }
   return(list(data=data_transformed,model=model_out))
}
# 2007-12-10:
#
#
compute_pattern <- function(data,fun=mean){
	data <- data.matrix(aggregate(data,fun,by=list(class=data$class),na.rm=TRUE))[,-1]
	row.names(data) <- 1:nrow(data)	
	if(!is.na(match("class", dimnames(data)[[2]])))
		data <- data[,-match("class", dimnames(data)[[2]])]
	return(data)
}
# 2008-03-22:
#
#
retrieve_model <-  function(canalysis,query){
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
# 2008-03-27:
#
#
fullproba_ftable <- function(z1, z2){
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
# 2008-03-22:
#
#
cross_compare_models <- function(model1,model2,fileNameCSV=today()){
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
# 2007-11-xx:
#
#
compare_canalysis <- function(analysis1,analysis2=NULL,query=NULL,filePrefix=NULL,bbox_threshold = 5){
   canalysis1 <- analysis1[["cluster_analysis"]]
   cross_comparison_results <- list()
   model      <- list()
   if(!is.null(analysis2) && !is.null(query)){
      # WE COMPARE THE RESULT OF TWO ANALYSIS FOR A GIVEN MODEL
      canalysis2 <- analysis2[["cluster_analysis"]]
      model[[1]] <- retrieve_model(canalysis1,query)
      model[[2]] <- retrieve_model(canalysis2,query)
      if(!is.null(filePrefix))
         filePrefix <- paste2(filePrefix,"_X_comparison_",gsub(",","_",query),".csv")
      if(!is.na(model[[1]][["model"]]$loglik) && !is.na(model[[2]][["model"]]$loglik))
         cross_comparison_results[[1]] <- cross_compare_models(model[[1]],model[[2]],fileNameCSV=filePrefix)
   }
   else{
      # SELECT ALL 2-BY-2 COMPARISONS BETWEEN MODELS THAT ARE WITHIN THE RELATIVE
      # 5% BBOX TO THE BEST BIC VALUE
      vertices <- get_vertice_set(canalysis1[["bic_table_relative"]], bbox_threshold)
      if(nrow(vertices)>=1){
         # PROCEED TO 2-BY-2 COMPARISON OVER THE SET OF VERTICES
         for(v in 1:nrow(vertices)){
            # RETRIEVE CLUSTERING Z-MAPS FROM 'CANALYSIS'
            model[[1]] <- retrieve_model(canalysis1,vertices[v,1])
            model[[2]] <- retrieve_model(canalysis1,vertices[v,2])
            if(!is.null(filePrefix))
               filePrefix <- paste2(filePrefix,"_",gsub(",","_",vertices[v,1]),"-",gsub(",","_",vertices[v,2]),".csv")
            if( !is.na(model[[1]][["model"]]$loglik) && !is.na(model[[2]][["model"]]$loglik))
               cross_comparison_results[[v]] <- cross_compare_models(model[[1]],model[[2]],fileNameCSV=filePrefix)
         }
      }
   }
   return(cross_comparison_results)
}
# 2008-01-04: model_based_clustering
#
#
model_based_clustering <- function(data,filePrefix=filePrefix,stats_fun,G_set=1:9
   #, clust_fun = em_clustering  
   , clust_fun = mclust2
   , model_names=c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV")
#        model_names=c("EII","VEI")
#    ,    model_names=c("EII","VVI")

        ) {
   # SUBSELECT OUTCOMES ON WHICH WE PROCEED TO THE CLUSTER ANALYSIS
   local_data                      <- data[["transformed"]][,data[["analysis_outcomes"]]]
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
            additional_factors  <- which( !(colnames(data[["transformed"]]) %in% colnames(local_data)))
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
                              canalysis_outcomes = data[["analysis_outcomes"]])
            for(i in names(stats_fun)){
               mbc[["out"]][[l]][["stats"]][[i]] <- stats_fun[[i]](data=data,class=mbc[["out"]][[l]][["z"]])
               cat(mbc[["out"]][[l]][["G"]],mbc[["out"]][[l]][["modelName"]],", ")
            }
         }
      }
   }
   return(mbc)
}
# 2008-04-06:
# START MODEL BASED CLUSTERING WITH AN EM STEP
#
em_mbc          <- function(data, modelName="VVI", G=5) {
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
# 2008-04-06:
# START MODEL BASED CLUSTERING WITH AN EM STEP
#
em_clustering          <- function(data, G_set=1:9,
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
# 2007-12-05: 
#
#
cluster_analysis <- function(data,filePrefix=filePrefix,G_set=1:9,
	stats_fun = list(       oddratios=function(data,class){return(statistics_logodds(data,class))}
#				,gen_naive_bayes=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_naive_bayes,predict=predict_naive_bayes),K=5))}
#				,gen_knn=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_knn,predict=predict_knn),K=5))}
#				,gen_svm=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_svm,predict=predict_svm),K=5))}
				)){
	# 
	# PROCEED TO MBC ANALYSIS ON THE TRANSFORMED DATA
	# 
	clustering <- model_based_clustering(data,filePrefix=filePrefix,stats_fun=stats_fun,G_set=G_set)	
	return(clustering)
}
# 2008-01-07:
#
#
statistical_patterns <- function(data,fun_list=list(avg=mean,median=median)){
        data_out <- list()
        for(fun in names(fun_list))
                data_out[[fun]] <- compute_pattern(data,fun=fun_list[[fun]])
        return(data_out)
}
# 2007-11-25: PROCEDURE DOING BOTH THE COMPUTATION OF THE 
# MODEL BASED CLUSTER ANALYSIS, OF THE AVG PATTERNS, AND 
# THE DISPLAY
# source('./SubtypeDiscovery.R') ; debug(analysis) ; debug(graphic_characterization); analysis(data=m[complete.cases(m),],analysis_outcomes = outcomes.parkinson,fun_transform=list(duration_random_effect=transform_remove_effect,min=transform_MIN,max=transform_MAX),filePrefix="GARP") ; 
analysis <- function(data,
	analysis_outcomes = list(),     # OUTCOMES ON WHICH THE TRANSFORMATION AND CLUSTERING SHOULD BE DONE
	sumscore_groups = list(),       # LOGODDS AND EVENTUALLY OTHER STATISTICS USE THIS GROUPING TO DO SUM SCORES
        stats_fun = list(),             # A LIST OF FUNCTION TO EVALUATE THE CLUSTERING RESULT
        G_set=1:9,
        formula_matrix = NULL,
	fun_transform=list(identity=function(x){return(data=I(x),model=NULL)}),filePrefix=NULL,...){
	#
	# EVENTUALLY INITIALIZE THE FIRST PART OF THE FILE-NAMING 
	# AND THE LIST OF OUTCOMES ON WHICH THE ANALYSIS WILL BE PERFORMED
	#
	if(!is.null(filePrefix)){
            filePrefix <- paste2(today(),"_",filePrefix)
            postscript(file=paste2(filePrefix,"Cluster_Analysis.ps"))
            }
	if(is.null(analysis_outcomes))
            analysis_outcomes <- dimnames(data)[[2]]
	#
	# 0. DATA SELECTION
	#		CONFINE TO COMPLETE CASES AND PROCEED TO THE ANALYSIS
	#
	data <- list(   orig    = data,
                        cc      = data[row.names(na.omit(data[,analysis_outcomes])),],
                        analysis_outcomes  = sort(analysis_outcomes),
                        sumscore_groups    = sumscore_groups,
                        formula_matrix     = formula_matrix,
                        mask_canalysis     = colnames(formula_matrix) %in% analysis_outcomes,
                        mask_evaluation    = colnames(formula_matrix) %in% unlist(sumscore_groups))
	#
	# 0. DATA SELECTION
	#		PROCEED TO AN ANALYSIS ON THE OUTCOMES 
	#	. CORRELATION
	#	. BOXPLOTS
	#	. HISTOGRAMS
	# 
        orig_data <- list(data=data[["cc"]],model=NULL,canalysis_outcomes=analysis_outcomes)
        outcomes_analysis(orig_data,filePrefix=paste2(filePrefix,"Original_Data_CC"))
	#
	# 1. AND 2. MULTIPLE CLUSTERING ANALYSIS AND MODEL SELECTION (BIC)
	#	. MODEL BASED CLUSTERING FROM 1-9 GROUPS WITH DIFFERENT COVARIANCE MATRICES
	#	. BIC TABLE FOR MODEL AND NUMBER OF GROUP SELECTION
	#	. DATA TRANSFORMATION
	# 
        transform_out <- transform_outcomes(data,fun_transform=fun_transform)
	data[["transformed"]] <- transform_out[["data"]]
	data[["fun_transform"]] <- fun_transform
	#
	# ALSO ANALYSE OUTCOMES FOR TRANSFORMED DATA
	#
	# outcomes_analysis(data[["transformed"]][,analysis_outcomes],filePrefix=paste2(filePrefix,"Transformed_Data"))
        transform_out[["canalysis_outcomes"]] <- analysis_outcomes
        outcomes_analysis(transform_out,filePrefix=paste2(filePrefix,"Transformed_Data"))
	#
	canalysis               <- cluster_analysis(data,stats_fun = stats_fun, filePrefix=filePrefix,G_set=G_set)
	data[["cc"]]            <- cbind(data[["cc"]],class=map(canalysis$z))
	data[["transformed"]]   <- cbind(data[["transformed"]],class=map(canalysis$z))
	#
	# CHARACTERIZE EACH CLUSTER BY ITS AVG PATTERN 
        #
	#
	# 3. ASSESS THE CLUSTERING METHOD CONSISTENCY
	#	DO MODEL COMPARISON ROW- AND COLUMN-WISE IN THE 5% BOUNDING BOX 
	# 
        data[["cluster_analysis"]] <- canalysis
	compare_canalysis(data,filePrefix=filePrefix,bbox_threshold=2)
        save_models(data,filePrefix=filePrefix)
	#
	# 4. CHARACTERIZE BY VISUALIZATION THE DIFFERENT CLUSTERS
	#	. PARALLEL COORDINATES (USER-DEFINED)
	#	. HEATMAPS (DEFAULT)
	#	. DENDROGRAMS ILLUSTRATING AVG-PATTERNS AND OUTCOME'S SIMILARITY
	# 
	if(!is.null(filePrefix))
                graphics.off()
	return(data)
}
#
# 2008-03-31
#
stability_assessment <- function(
     data               = na.omit(as.data.frame(df))
   , analysis_outcomes  = NULL
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
                                                transform_addnoise(data, model, canalysis_outcomes = analysis_outcomes, type, 
                                                relative_noise_amount = noise_vector[noise_idx], rand_seed = 6013+3*i))}
         filePrefix_local <- paste2(filePrefix,noise_vector[noise_idx],"_",adj,"_")
         # PROCEED TO CLUSTER ANALYSIS 10 TIMES WITH A DIFFERENT SEED
         for(i in 1:10)
            canalysis_results[[adj]][[noise_idx]][[i]] <- analysis(
                                                              data                   = data
                                                            , analysis_outcomes      = analysis_outcomes 
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
###################### END, FUNCTION DEFINITION SECTION ####################
# 2008-01-13: 
#
#
get_model_set <- function(relativeBic, bbox_threshold=5){
   relativeBic <- list(bicTable=relativeBic, bbox = apply(abs(relativeBic) < bbox_threshold,2,which))
   relativeBic$bbox <- relativeBic$bbox[which(lapply(relativeBic$bbox,length) >0)]
   # DERIVE VERTICES TO COMPARE WITHIN THE 5% BBOX OF THE BIC TABLE
   relativeBic$best_models <- matrix(0,0,2,dimnames=list(list(),list("modelName","G")))
   for(m in names(relativeBic$bbox))
      for(g in relativeBic$bbox[[m]])
         relativeBic$best_models <- rbind(relativeBic$best_models,c(m,row.names(relativeBic[["bicTable"]])[g]))
   return(relativeBic[["best_models"]])
}
# 2007-12-03: 
# FROM A TABLE OF RELATIVE BIC VALUES, RETRIEVE ALL THE POSSIBLE PAIRS OF MODELS
# TO DO COMPARISON ON. THE MODELS ARE LESS THAN (%) BBOX_THRESHOLD 'WORST' THAN THE OPTIMAL 
# BIC VALUE. 
get_vertice_set <- function(relativeBic, bbox_threshold=5){
   relativeBic <- list( bicTable=relativeBic, 
                        # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
                        bbox = apply(relativeBic <= max(sort(relativeBic)[1:bbox_threshold]), 2,which)
                        # TO SELECT ALL MODELS BELOW A GIVEN THRESHOLD
                        # bbox = apply(abs(relativeBic) < bbox_threshold,2,which)
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
###########################################################
# STATISTICAL EVALUATION OF THE PATTERNS PROVIDED THE     #
# POPULATION DISTRIBUTION                                 #
###########################################################
# 2006-11-23:
# THIS FUNCTION DERIVE LAMBDA SIBS FROM A SIB PAIR FAMILY DATASET
# OTHER STATISTICS ARE DERIVED SUCH AS:
# - POPULATION WIDE CLUSTER PREVALENCE
# - PROBANT / SIBLING CLUSTER PREVALENCE
# - SIB PAIRS CONJONCTIVE PROBABILITIES
# - AND LAMBDA SIBS
# get as input a matrix with the "family" as rows and sib1/probant, sib2 as columns
statistics_lambdasibs <- function(data,class)
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
# 2007-12-04:
#
#
statistics_logodds <- function(data,class,fun_midthreshold=median){
        #
        class           <- map(class)
	local_data      <- cbind(data[["transformed"]],class=class)
	s               <- as.data.frame(matrix(0,length(unique(local_data$class)),length(data[["sumscore_groups"]])))
	# FOR EACH HIERARCHICAL SUBSET OF OUTCOMES, MAKE SUM SCORES
	for(gr in names(data[["sumscore_groups"]])){
		mat.l   <- cbind(SScore=apply(data.matrix(local_data[,data[["sumscore_groups"]][[gr]]]),1,sum),class=local_data$class)
		# FOR EACH PATTERN, DERIVE ITS STATISTIC
		for(g in sort(unique(local_data$class))){
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
			s[g,match(gr,names(data[["sumscore_groups"]]))] <- log(mat[1,1] * mat[2,2] / (mat[1,2] * mat[2,1]))
			dimnames(s)[[2]][match(gr,names(data[["sumscore_groups"]]))] <- paste(gr,sprintf("%.1f",med_SScore),sep="_",collapse="")
		}
	}
	return(list(out=s))
}
###########################################################
# SAMPLING FUNCTION TO PREPARE TRAIN/TEST SPLIT           #
###########################################################
# 2008-01-10
#
#
stratified_traintest_split <- function(data,K=10,perc=0.7){
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
# 2007-12-04:
# EXTENDS FUNCTION GENERATE.CV FROM WILCOXCV PACKAGE
#
stratified_cv_folds <- function(n, K, G) 
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
# 2007-12-10
#
#
stratified_cv <- function(data,K=10){
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
###########################################################
# REPRODUCIBILITY: AI INTERFACES AND EVALUATION           #
###########################################################
# 2007-12-10
#
#
model_naive_bayes <- function(formula, trainset, testset, laplace = 1){
   return(naiveBayes(class ~ . , data=trainset, laplace = 1))
}
# 2007-12-10
#
#
model_knn <- function(formula, trainset, testset, k = 1, l = 0, prob = FALSE, use.all = TRUE){
   return(knn(trainset[,-match("class",names(trainset))], testset[,-match("class",names(testset))] ,trainset[,"class"], k = k, l = l, prob = prob, use.all = use.all))
}
# 2007-12-10
#
#
model_svm <- function(formula, trainset, testset, kernel="linear",type="C-classification"){
   return(svm(class ~ . , data=trainset,kernel=kernel,type=type))
} 
# 2007-12-10
#
#
predict_naive_bayes <- function(model,testset){
   return(map(predict(model,testset[,-match("class",names(testset))],type="raw")))
}
# 2007-12-10
#
#
predict_knn <- function(model, testset){
   return(model)
}
# 2007-12-10
#
#
predict_svm <- function(model, testset){
   predict(model,testset[,-match("class",names(testset))])
}
# 2007-12-04: 
# 
# 
evaluate_generalization <- function(data,class
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
            trainset        <- cbind(fun_out[["data"]][,data[["analysis_outcomes"]]],class=trainset$class)
            fun_out         <- fun_transform[[fun_name]](data=testset[,-match("class",dimnames(testset)[[2]])],model=fun_out[["model"]])
            testset         <- cbind(fun_out[["data"]][,data[["analysis_outcomes"]]],class=testset$class)
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
###########################################################
# OTHER STATS                                             #
###########################################################
area_under_uncertainty_curve <- function(class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}
# 
#
#
missing_at_random <- function(data,class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}
# 2008-01-13: 
#
#
get_model_set <- function(relativeBic, bbox_threshold=5){
   relativeBic <- list(bicTable=relativeBic, bbox = apply(abs(relativeBic) < bbox_threshold,2,which))
   relativeBic$bbox <- relativeBic$bbox[which(lapply(relativeBic$bbox,length) >0)]
   # DERIVE VERTICES TO COMPARE WITHIN THE 5% BBOX OF THE BIC TABLE
   relativeBic$best_models <- matrix(0,0,2,dimnames=list(list(),list("modelName","G")))
   for(m in names(relativeBic$bbox))
      for(g in relativeBic$bbox[[m]])
         relativeBic$best_models <- rbind(relativeBic$best_models,c(m,row.names(relativeBic[["bicTable"]])[g]))
   return(relativeBic[["best_models"]])
}
# 2007-12-03: 
# FROM A TABLE OF RELATIVE BIC VALUES, RETRIEVE ALL THE POSSIBLE PAIRS OF MODELS
# TO DO COMPARISON ON. THE MODELS ARE LESS THAN (%) BBOX_THRESHOLD 'WORST' THAN THE OPTIMAL 
# BIC VALUE. 
get_vertice_set <- function(relativeBic, bbox_threshold=5){
   relativeBic <- list( bicTable=relativeBic, 
                        # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
                        bbox = apply(relativeBic <= max(sort(relativeBic)[1:bbox_threshold]), 2,which)
                        # TO SELECT ALL MODELS BELOW A GIVEN THRESHOLD
                        # bbox = apply(abs(relativeBic) < bbox_threshold,2,which)
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
###########################################################
# STATISTICAL EVALUATION OF THE PATTERNS PROVIDED THE     #
# POPULATION DISTRIBUTION                                 #
###########################################################
# 2006-11-23:
# THIS FUNCTION DERIVE LAMBDA SIBS FROM A SIB PAIR FAMILY DATASET
# OTHER STATISTICS ARE DERIVED SUCH AS:
# - POPULATION WIDE CLUSTER PREVALENCE
# - PROBANT / SIBLING CLUSTER PREVALENCE
# - SIB PAIRS CONJONCTIVE PROBABILITIES
# - AND LAMBDA SIBS
# get as input a matrix with the "family" as rows and sib1/probant, sib2 as columns
statistics_lambdasibs <- function(data,class)
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
# 2007-12-04:
#
#
statistics_logodds <- function(data,class,fun_midthreshold=median){
        #
        class           <- map(class)
	local_data      <- cbind(data[["transformed"]],class=class)
	s               <- as.data.frame(matrix(0,length(unique(local_data$class)),length(data[["sumscore_groups"]])))
	# FOR EACH HIERARCHICAL SUBSET OF OUTCOMES, MAKE SUM SCORES
	for(gr in names(data[["sumscore_groups"]])){
		mat.l   <- cbind(SScore=apply(data.matrix(local_data[,data[["sumscore_groups"]][[gr]]]),1,sum),class=local_data$class)
		# FOR EACH PATTERN, DERIVE ITS STATISTIC
		for(g in sort(unique(local_data$class))){
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
			s[g,match(gr,names(data[["sumscore_groups"]]))] <- log(mat[1,1] * mat[2,2] / (mat[1,2] * mat[2,1]))
			dimnames(s)[[2]][match(gr,names(data[["sumscore_groups"]]))] <- paste(gr,sprintf("%.1f",med_SScore),sep="_",collapse="")
		}
	}
	return(list(out=s))
}
###########################################################
# SAMPLING FUNCTION TO PREPARE TRAIN/TEST SPLIT           #
###########################################################
# 2008-01-10
#
#
stratified_traintest_split <- function(data,K=10,perc=0.7){
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
# 2007-12-04:
# EXTENDS FUNCTION GENERATE.CV FROM WILCOXCV PACKAGE
#
stratified_cv_folds <- function(n, K, G) 
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
# 2007-12-10
#
#
stratified_cv <- function(data,K=10){
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
###########################################################
# REPRODUCIBILITY: AI INTERFACES AND EVALUATION           #
###########################################################
# 2007-12-10
#
#
model_naive_bayes <- function(formula, trainset, testset, laplace = 1){
   return(naiveBayes(class ~ . , data=trainset, laplace = 1))
}
# 2007-12-10
#
#
model_knn <- function(formula, trainset, testset, k = 1, l = 0, prob = FALSE, use.all = TRUE){
   return(knn(trainset[,-match("class",names(trainset))], testset[,-match("class",names(testset))] ,trainset[,"class"], k = k, l = l, prob = prob, use.all = use.all))
}
# 2007-12-10
#
#
model_svm <- function(formula, trainset, testset, kernel="linear",type="C-classification"){
   return(svm(class ~ . , data=trainset,kernel=kernel,type=type))
} 
# 2007-12-10
#
#
predict_naive_bayes <- function(model,testset){
   return(map(predict(model,testset[,-match("class",names(testset))],type="raw")))
}
# 2007-12-10
#
#
predict_knn <- function(model, testset){
   return(model)
}
# 2007-12-10
#
#
predict_svm <- function(model, testset){
   predict(model,testset[,-match("class",names(testset))])
}
# 2007-12-04: 
# 
# 
evaluate_generalization <- function(data,class
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
            trainset        <- cbind(fun_out[["data"]][,data[["analysis_outcomes"]]],class=trainset$class)
            fun_out         <- fun_transform[[fun_name]](data=testset[,-match("class",dimnames(testset)[[2]])],model=fun_out[["model"]])
            testset         <- cbind(fun_out[["data"]][,data[["analysis_outcomes"]]],class=testset$class)
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
###########################################################
# OTHER STATS                                             #
###########################################################
area_under_uncertainty_curve <- function(class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}
# 
#
#
missing_at_random <- function(data,class){
   # AREA UNDER THE UN-CERTAINTY CURVE FOR EACH TWO MODELS
   auuc         <- mean(1-apply(class,1,max,na.rm=TRUE))
   auuc_sd      <- sd(  1-apply(class,1,max,na.rm=TRUE))
   out          <- matrix(NA,ncol(class),2,dimnames=list(1:ncol(class),c("AUUC","AUUC sdev")))
   out[1,1]     <- auuc
   out[1,2]     <- auuc_sd
   return(list(auuc = auuc, sd = auuc_sd, out = out))
}
# 2007-11-24: DISPLAY OF MULTIPLE HISTOGRAMS ON EACH PAGE
# FOR COMPARISON PURPOSE WHEN APPLYING SOME TRANSFORMATION
# ON THE OUTCOMES
graphic_histograms <- function(cdata,fileName="2007-11-24_IMG_distributions.ps",xlim=NULL){
        local_data <- cdata[["data"]][,cdata[["canalysis_outcomes"]]]
        local_data <- as.data.frame(local_data)
	par('mai'=c(0.5,0.5,0.5,1),'mfrow'=c(4,5))
        #
        hist_list  <- list()
        hist_stats <- list()
        hist_counts <- c()
        models <- cdata[["model"]]
        for(i in colnames(local_data)){
           hist_list[[i]]  <- hist(local_data[,i],plot=FALSE)
           local_stat_txt <- ""
           hist_stats[[i]] <- list(txt="",xlim=list(),ylim=list())
           for(j in names(models)){
               local_model <- models[[j]][["model"]][[i]]
               if(models[[j]][["type"]] == "lm")
                  local_stat_txt <- paste2(local_stat_txt,local_model[["formula"]],"\n")
               else
                  local_stat_txt <- paste2(local_stat_txt,sprintf("%1.2f",local_model[[1]]),", ")
           }
           local_stat_txt <- paste2(local_stat_txt,length(table(local_data[,i]))," distinct values")
           hist_stats[[i]][["txt"]] <- local_stat_txt
           hist_stats[[i]][["xlim"]] <- c(min(local_data[,i],na.rm=TRUE),max(local_data[,i],na.rm=TRUE))
#           hist_stats[[i]][["xlim"]] <- c(quantile(data.matrix(local_data[,i]),probs=c(0.05),na.rm=TRUE),
#                quantile(data.matrix(local_data[,i]),probs=c(0.95),na.rm=TRUE))
           hist_stats[[i]][["ylim"]] <- range(hist_list[[i]][["counts"]])
        }
        #
	for(i in 1:ncol(local_data)){
            par('lab'=c(3,3,3),  # number of ticks
                'mar'=c(2.5,2.5,1.5,1) # margins bottom,left,top,right
                )
            plot(hist_list[[i]],
                main=dimnames(local_data)[[2]][i],
                xlab="",ylab="",xlim=hist_stats[[i]][["xlim"]],ylim=hist_stats[[i]][["ylim"]]
                #,ylog=TRUE
                ) #!!!!
            text(x=mean(hist_stats[[i]][["xlim"]]),y=hist_stats[[i]][["ylim"]][2]*0.85,hist_stats[[i]][["txt"]])
            }
		
}
# 2007-12-05:
#
#
outcomes_analysis <- function(cdata,model=NULL,canalysis_outcomes,filePrefix=today(),xlim=c(-3,3)){
        local_data <- cdata[["data"]]
	#
	# CORRELATION
	#	. CORRELATION MATRIX
	#	. ABSOLUTE VALUED CORRELATION MATRIX TAKEN AS A DISTANCE MATRIX
	#	TO PLOT A DENDROGRAM WHICH ILLUSTRATES OUTCOME'S RELATION
	#
	cor_mat <- cor(local_data[,cdata[["canalysis_outcomes"]]],use="pairwise.complete.obs")
	write.table(cor_mat,file=paste2(filePrefix,"Correlation_Matrix.csv"),dec=",",sep=";")
	hclust(as.dist(abs(cor_mat)))
	#
	# BOXPLOTS 
	#		SUMMARY DITRIBUTIONAL SHAPE OF THE OUTCOMES BY THEIR BOX PLOTS
	#	. 1ST AND 3RD QUARTILES
	#	. MEDIAN
	#	. 95% BOUNDS
	#	. + EXTREMES BEYOND THE 95%
	# 
	par(mfrow=c(1,1),'mai'=c(1.2,0.7,0.7,0.7),new=FALSE)
	boxplot(local_data,'las'=2,main=paste2(filePrefix,"Boxplot"))
	#
	# HISTOGRAMS 
	#
	graphic_histograms(cdata,model,xlim=xlim)
}
# 2006-12-15: (ORIGINALY)
# 2007-03-16: ADAPTED TO TAKE INTO ACCOUNT BOOSTRAP 95%CI
# INPUT: GET AS INPUT AN ARRAY OF 
# 		(OBSERVED MEAN PATTERN: COLUMN VECTORS BEING THE CLUSTER CENTROIDS)
#		(BOOTSTRAP MEAN PATTERN: COLUMN VECTORS OF BOOTSTRAP CLUSTER CENTROID)
#		(95% LOWER BOUND CI: COLUMN VECTOR ESTIMATED BY BOOTSTRAP ON EACH SITE )
# 		(95% UPPER BOUND CI: COLUMN VECTOR ESTIMATED BY BOOTSTRAP ON EACH SITE )
# 		(FOR EACH CLUSTER, FOR EACH MODEL)
# OUTPUT: 	VIGNETTES
#			SUMMARY STATISTIC OF THE SUM OF 95% COVERAGES FOR EACH MODEL
# DO:
graphic_characterization	<- function(data, class_count=class_count,
	graphics=graphic_params,
	titlePrefix=NULL,
        canalysis_outcomes = canalysis_outcomes){
        #
        local_data <- data 
	# LOAD TRANSFORMATION MATRIX (X,Y,Z)
	visu_mat <- read.csv(file_longitudinal_visu,header=TRUE,row.names=1,dec=",")
	# visu_mat <- visu_mat[complete.cases(visu_mat[dimnames(data[[1]])[[2]],]),]
	vignettes.title	<- levels(visu_mat$Z)[levels(visu_mat$Z) != 0]
        # WATCH-OUT! THE PATTERN [[1]] IS EXPECTED TO BE THE ONE FROM WHICH THE
        # ORDERING IS DERIVED THUS, THESE SHOULD REFER TO THE MEAN OR THE
        # MEDIAN PATTERNS
        color_sel       <- get_coloring_scheme(data[[1]],canalysis_outcomes = canalysis_outcomes)
	par(mfrow=c(3,3), 'las'=1,'mai'=c(0.4,0.55,0.18,0.05))
	#
	# LOOP OVER EACH VIGNETTE
	#
	for(v in vignettes.title){
		# TAKE THE SUBSET OF PHENOTYPES TO BE VISUALIZED ON VIGNETTE Z: TRANSFORMATION MATRIX
		T_s					<- visu_mat[visu_mat$Z == v,]
		plot(x=0        ,new=TRUE
                                ,ann=FALSE
                                ,pch=18
                                ,xlim=c(graphics[["xlim"]]$min,graphics[["xlim"]]$max)
                                ,ylim=c(graphics[["ylim"]]$min,graphics[["ylim"]]$max)
                                ,col="white",axes=FALSE)
                title(main=paste2(titlePrefix,v))
		# LOOP OVER EACH GROUP PATTERN THAT IS TO VISUALIZE
		for(g in color_sel[["cluster_order"]]){
			# LOOP OVER THE DIFFERENT STATISTICS FOR EACH PATTERN
			for(s in names(data)){
                                local_pattern           <- data[[s]][color_sel[["cluster_order"]],]
				C_i_s                   <- cbind(as.numeric(local_pattern[g,dimnames(T_s)[[1]]]),0)
				D_i_s			<- C_i_s+T_s[,1:2]
				D_i_s			<- D_i_s[sort.list(D_i_s$Y),]
                                gap                     <- D_i_s[2,"Y"] - D_i_s[1,"Y"]
				for(l in 1:(nrow(D_i_s)-1)){
                                        is_white_gap <- (D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap)
					if((s == "median" || s == "avg") & !is_white_gap )
						arrows(D_i_s[l,1],D_i_s[l,2],D_i_s[l+1,1],D_i_s[l+1,2],col=color_sel[["cluster_color"]][g],length=0,lwd=3)
                                        if((s == "lowquant" || s == "upquant") & !is_white_gap)
						arrows(D_i_s[l,1],D_i_s[l,2],D_i_s[l+1,1],D_i_s[l+1,2],col=color_sel[["cluster_color"]][g],length=0,lty="dashed",lwd=1)
                                        # else, as (is_white_gap == TRUE) then do not draw any arrow...
				}
				if(v == graphics$legend_on){
					# PLOT THE LEGEND ON ONE OF THE VIGNETTES			
					# DO THE LEGEND WITHIN [3.5;6.5] AND START FROM THE TOP, I.E. 6.5
					#y1 <- y0 <- 0.9*graphics[["ylim"]]$max-(g-1)*3/nrow(data[[1]])
					y1 <- y0 <- graphics[["legend_ycoord"]] + 5*g/nrow(data[[1]])
					x0 <- graphics[["legend_xcoord"]]
					x1 <- 1.2*graphics[["legend_xcoord"]]
					arrows(x0,y0,x1,y1,col=color_sel[["cluster_color"]][g],length=0,lwd=3)
                                        local_g <- color_sel[["cluster_order"]][g]
					text(x0+0.15,y1,labels=paste2(local_g," (",class_count[as.numeric(local_g)],")"),pos=4)
				}						
			}
			# DIMNAMES
			axis(2,at=D_i_s[,2],labels=dimnames(D_i_s)[[1]],las=2,tick=FALSE)
			axis(1,at=seq(from = graphics[["xlim"]]$min , to = graphics[["xlim"]]$max , by = ((graphics[["xlim"]]$max - graphics[["xlim"]]$min)/4)))
		}
	}
	# DEFINE COLOR GRADIENT
	gradient.colors <- graphics[["color_gradient"]] 
	# DETERMINES ROW AND COLUMN ORDERINGS OF THE AVG PATTERN
        data <- data[[1]][,canalysis_outcomes]
	rowD <- hclust(dist(data))
	colD <- hclust(dist(t(data)))
	# d <- data[rowD$order,colD$order]
        d <- data[rowD$order, rev(sort(colnames(data)))]
	# PRODUCE THE 3-PLOTS: IMAGE/HEATMAP, ROW AND COLUMN DENDROGRAMS
	image(1L:nrow(d), 1L:ncol(d), d, xlim = 0.5 + c(0, nrow(d)), ylim = 0.5 + c(0, ncol(d)), axes = FALSE, xlab = "", ylab = "",col=gradient.colors,main="Average pattern visualization")
	axis(2, 1L:ncol(d), labels = dimnames(d)[[2]], las = 2, line = -0.5, tick = 0, cex.axis = graphics[["cex"]])
	axis(1, 1L:nrow(d), labels = dimnames(d)[[1]], las = 1, line = -0.5, tick = 0, cex.axis = graphics[["cex"]])
	plot(rowD, axes = FALSE, yaxs = "i", main="Patterns arranged by similarity",ylab=NULL,cex=graphics[["cex"]])
	plot(colD, axes = FALSE, yaxs = "i", main="Outcomes arranged by similarity",ylab=NULL,cex=graphics[["cex"]])
}
# 2008-02-18: Exploratory Data Analysis for longitudinal data, profile plots
#
#
graphic_profile <- function(data,psoutput=NULL,plot_sample=20){
   if(!is.null(psoutput))
      postscript(file=psoutput,paper="a4")
   par(mfrow=c(3,3), 'las'=1,'mai'=c(0.3,0.3,0.3,0.5))
   id_set    <- sample(row.names(data))
   #
   for(o in colnames(data[,,1])){
      a <- as.data.frame(data[,o,])
      b <- reshape(a,direction="long",ids=row.names(a),times=1:dim(data)[[3]],timevar="Year",varying=list(names(a)),v.names=o)
      # INIT PLOT
      plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(min(a,na.rm=TRUE),max(a,na.rm=TRUE)),ylab=o)
      # RAW DATA
      for(i in 1:length(id_set)){
         if((length(id_set)-i) < plot_sample)
            line_col <- "blue"
         else
            line_col <- "black"
         lines(xy.coords(b[b$id == id_set[i],c("Year",o)]),type="l",col=line_col)
         }
      # STATS
      for(s in list(list(fun=function(x){return(mean(x,na.rm=TRUE))},col="red",lty="solid",lwd = 2),
                  list(fun=function(x){return(median(x,na.rm=TRUE))},col="red",lty="dashed",lwd = 1),
                  list(fun=function(x){return(quantile(x, na.rm=TRUE, probs = seq(0, 1, 0.05))["5%"])},col="red",lty="dashed",lwd = 1),
                  list(fun=function(x){return(quantile(x, na.rm=TRUE, probs = seq(0, 1, 0.05))["95%"])},col="red",lty="dashed",lwd = 1)
                  )){
         s[["pattern"]] <- cbind(Year=1:dim(data)[[3]],o=apply(data[,o,],2,s[["fun"]]))
         lines(xy.coords(s[["pattern"]]),type="l",col = s[["col"]],lty = s[["lty"]],lwd = s[["lwd"]])
      }
   }
   if(!is.null(psoutput))
      graphics.off()
}
# 2008-02-18: Exploratory Data Analysis for longitudinal data, profile plots
#
#
graphic_cluster_profiles <- function(data,canalysis_outcomes = parkinson_canalysis_outcomes,
        psoutput        =NULL,
        fun_stats       =list(
#                f_median   = list(fun = function(x){return(median(x,na.rm=TRUE))},lty="solid",lwd = 4)
                 f_mean   = list(fun = function(x){return(mean(x,na.rm=TRUE))},  lty="solid",lwd = 4)
                , f_25quant = list(fun = function(x){return(quantile(x,probs=0.25,na.rm=TRUE))},  lty="dashed",lwd = 1)
                , f_75quant = list(fun = function(x){return(quantile(x,probs=0.75,na.rm=TRUE))},  lty="dashed",lwd = 1)
                ))
                {
   #
   if(!is.null(psoutput))
      postscript(file=psoutput,paper="a4")
   #
   par(mfrow=c(2,9), 'las'=1,'mai'=c(0.3,0.1,0.3,0.3),oma=c(1,1,1,2))
   avg_pattern  <- compute_pattern(as.data.frame(data[,c(canalysis_outcomes,"class"),1]),fun = fun_stats[[1]][[1]])
   color_sel    <- get_coloring_scheme(avg_pattern,canalysis_outcomes = canalysis_outcomes)
   ymin         <- min(quantile(data,probs=0.05,na.rm=TRUE),na.rm=TRUE)
   ymax         <- max(quantile(data,probs=0.95,na.rm=TRUE),na.rm=TRUE)
   #
   id_set    <- sample(row.names(data))
   for(o in  canalysis_outcomes){
      # INIT PLOT
      plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(ymin,ymax),main=o,ann=TRUE,axes=FALSE)
      # STATS
      for(cl in 1:length(color_sel[["cluster_order"]])){
         df_subset <- data[data[,"class",1] == color_sel[["cluster_order"]][cl],,]
         df_subset <- df_subset[row.names(df_subset)[!is.na(row.names(df_subset))],,]
         # RAW DATA
         a <- as.data.frame(df_subset[,o,])
         b <- reshape(a,direction="long",ids=row.names(a),times=1:dim(data)[[3]],timevar="Year",varying=list(names(a)),v.names=o)
         # RAW PROFILES
#         for(i in sample(unique(b$id))[1:5])
#            lines(xy.coords(b[b$id == i,c("Year",o)]),type="l",col=color_sel[["cluster_color"]][cl],lwd=1,lty="dotted")
         for(s in fun_stats){
            s[["pattern"]] <- cbind(Year=1:dim(data)[[3]],o=apply(df_subset[,o,],2,s[["fun"]]))
            lines(xy.coords(s[["pattern"]]),type="l",col = color_sel[["cluster_color"]][cl],lty = s[["lty"]],lwd = s[["lwd"]])
         }
         # PLOT THE LEGEND ON ONE OF THE VIGNETTES			
         # DO THE LEGEND WITHIN [3.5;6.5] AND START FROM THE TOP, I.E. 6.5
      }
      # Y-AXIS PLOT (2)
      axis(2,at=pretty(c(ymin,ymax),n=3))
   }
   #
   # LEGEND
   #
   plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(ymin,ymax),main="legend",ann=TRUE,axes=FALSE)
   for(cl in 1:length(color_sel[["cluster_order"]])){
      yrange <- (ymax-ymin)/length(color_sel[["cluster_order"]])
      y1 <- y0 <- ymax-cl*0.5*yrange
      x0 <- 1
      x1 <- 2
      arrows(x0,y0,x1,y1,col=color_sel[["cluster_color"]][cl],length=0,lwd=3)
      text(1.15*x1,y1,labels=paste2(cl," (",nrow(data[data[,"class",1] == color_sel[["cluster_order"]][cl],,]),")"),pos=4)
   }
   #
   #
   #
   if(!is.null(psoutput))
      graphics.off()
}
# 2008-02-25:
#
# TAKE AS INPUT A DATA ARRAY
get_coloring_scheme <- function(data, canalysis_outcomes = parkinson_canalysis_outcomes){
   data <- data.matrix(data)
   # MAKE A DENDROGRAM FROM THESE AVG PATTERN FOR THE DENDROGRAM ORDERING
   dendro_patterns      <- hclust(dist(data))
   # SELECT DIVERGING COLORS FROM RColorBrewer
   cluster_colors       <- brewer.pal(nrow(data),"Set1")
   # RETURN A DATA STRUCTURE WITH CLUSTERS RIGHTLY ORDERED AND THEIR ASSOCIATED COLOR
   return(list(cluster_order=dendro_patterns$order,cluster_color=cluster_colors))
}
# 2007-11-24: DISPLAY OF MULTIPLE HISTOGRAMS ON EACH PAGE
# FOR COMPARISON PURPOSE WHEN APPLYING SOME TRANSFORMATION
# ON THE OUTCOMES
graphic_histograms <- function(cdata,fileName="2007-11-24_IMG_distributions.ps",xlim=NULL){
        local_data <- cdata[["data"]][,cdata[["canalysis_outcomes"]]]
        local_data <- as.data.frame(local_data)
	par('mai'=c(0.5,0.5,0.5,1),'mfrow'=c(4,5))
        #
        hist_list  <- list()
        hist_stats <- list()
        hist_counts <- c()
        models <- cdata[["model"]]
        for(i in colnames(local_data)){
           hist_list[[i]]  <- hist(local_data[,i],plot=FALSE)
           local_stat_txt <- ""
           hist_stats[[i]] <- list(txt="",xlim=list(),ylim=list())
           for(j in names(models)){
               local_model <- models[[j]][["model"]][[i]]
               if(models[[j]][["type"]] == "lm")
                  local_stat_txt <- paste2(local_stat_txt,local_model[["formula"]],"\n")
               else
                  local_stat_txt <- paste2(local_stat_txt,sprintf("%1.2f",local_model[[1]]),", ")
           }
           local_stat_txt <- paste2(local_stat_txt,length(table(local_data[,i]))," distinct values")
           hist_stats[[i]][["txt"]] <- local_stat_txt
           hist_stats[[i]][["xlim"]] <- c(min(local_data[,i],na.rm=TRUE),max(local_data[,i],na.rm=TRUE))
#           hist_stats[[i]][["xlim"]] <- c(quantile(data.matrix(local_data[,i]),probs=c(0.05),na.rm=TRUE),
#                quantile(data.matrix(local_data[,i]),probs=c(0.95),na.rm=TRUE))
           hist_stats[[i]][["ylim"]] <- range(hist_list[[i]][["counts"]])
        }
        #
	for(i in 1:ncol(local_data)){
            par('lab'=c(3,3,3),  # number of ticks
                'mar'=c(2.5,2.5,1.5,1) # margins bottom,left,top,right
                )
            plot(hist_list[[i]],
                main=dimnames(local_data)[[2]][i],
                xlab="",ylab="",xlim=hist_stats[[i]][["xlim"]],ylim=hist_stats[[i]][["ylim"]]
                #,ylog=TRUE
                ) #!!!!
            text(x=mean(hist_stats[[i]][["xlim"]]),y=hist_stats[[i]][["ylim"]][2]*0.85,hist_stats[[i]][["txt"]])
            }
		
}
# 2007-12-05:
#
#
outcomes_analysis <- function(cdata,model=NULL,canalysis_outcomes,filePrefix=today,xlim=c(-3,3)){
        local_data <- cdata[["data"]]
	#
	# CORRELATION
	#	. CORRELATION MATRIX
	#	. ABSOLUTE VALUED CORRELATION MATRIX TAKEN AS A DISTANCE MATRIX
	#	TO PLOT A DENDROGRAM WHICH ILLUSTRATES OUTCOME'S RELATION
	#
	cor_mat <- cor(local_data[,cdata[["canalysis_outcomes"]]],use="pairwise.complete.obs")
	write.table(cor_mat,file=paste2(filePrefix,"Correlation_Matrix.csv"),dec=",",sep=";")
	hclust(as.dist(abs(cor_mat)))
	#
	# BOXPLOTS 
	#		SUMMARY DITRIBUTIONAL SHAPE OF THE OUTCOMES BY THEIR BOX PLOTS
	#	. 1ST AND 3RD QUARTILES
	#	. MEDIAN
	#	. 95% BOUNDS
	#	. + EXTREMES BEYOND THE 95%
	# 
	par(mfrow=c(1,1),'mai'=c(1.2,0.7,0.7,0.7),new=FALSE)
	boxplot(local_data,'las'=2,main=paste2(filePrefix,"Boxplot"))
	#
	# HISTOGRAMS 
	#
	graphic_histograms(cdata,model,xlim=xlim)
}
# 2006-12-15: (ORIGINALY)
# 2007-03-16: ADAPTED TO TAKE INTO ACCOUNT BOOSTRAP 95%CI
# INPUT: GET AS INPUT AN ARRAY OF 
# 		(OBSERVED MEAN PATTERN: COLUMN VECTORS BEING THE CLUSTER CENTROIDS)
#		(BOOTSTRAP MEAN PATTERN: COLUMN VECTORS OF BOOTSTRAP CLUSTER CENTROID)
#		(95% LOWER BOUND CI: COLUMN VECTOR ESTIMATED BY BOOTSTRAP ON EACH SITE )
# 		(95% UPPER BOUND CI: COLUMN VECTOR ESTIMATED BY BOOTSTRAP ON EACH SITE )
# 		(FOR EACH CLUSTER, FOR EACH MODEL)
# OUTPUT: 	VIGNETTES
#			SUMMARY STATISTIC OF THE SUM OF 95% COVERAGES FOR EACH MODEL
# DO:
graphic_characterization	<- function(data, class_count=class_count,
	graphics=graphic_params,
	titlePrefix=NULL,
        canalysis_outcomes = canalysis_outcomes){
        #
        local_data <- data 
	# LOAD TRANSFORMATION MATRIX (X,Y,Z)
	visu_mat <- read.csv(file_longitudinal_visu,header=TRUE,row.names=1,dec=",")
	# visu_mat <- visu_mat[complete.cases(visu_mat[dimnames(data[[1]])[[2]],]),]
	vignettes.title	<- levels(visu_mat$Z)[levels(visu_mat$Z) != 0]
        # WATCH-OUT! THE PATTERN [[1]] IS EXPECTED TO BE THE ONE FROM WHICH THE
        # ORDERING IS DERIVED THUS, THESE SHOULD REFER TO THE MEAN OR THE
        # MEDIAN PATTERNS
        color_sel       <- get_coloring_scheme(data[[1]],canalysis_outcomes = canalysis_outcomes)
	par(mfrow=c(3,3), 'las'=1,'mai'=c(0.4,0.55,0.18,0.05))
	#
	# LOOP OVER EACH VIGNETTE
	#
	for(v in vignettes.title){
		# TAKE THE SUBSET OF PHENOTYPES TO BE VISUALIZED ON VIGNETTE Z: TRANSFORMATION MATRIX
		T_s					<- visu_mat[visu_mat$Z == v,]
		plot(x=0        ,new=TRUE
                                ,ann=FALSE
                                ,pch=18
                                ,xlim=c(graphics[["xlim"]]$min,graphics[["xlim"]]$max)
                                ,ylim=c(graphics[["ylim"]]$min,graphics[["ylim"]]$max)
                                ,col="white",axes=FALSE)
                title(main=paste2(titlePrefix,v))
		# LOOP OVER EACH GROUP PATTERN THAT IS TO VISUALIZE
		for(g in color_sel[["cluster_order"]]){
			# LOOP OVER THE DIFFERENT STATISTICS FOR EACH PATTERN
			for(s in names(data)){
                                local_pattern           <- data[[s]][color_sel[["cluster_order"]],]
				C_i_s                   <- cbind(as.numeric(local_pattern[g,dimnames(T_s)[[1]]]),0)
				D_i_s			<- C_i_s+T_s[,1:2]
				D_i_s			<- D_i_s[sort.list(D_i_s$Y),]
                                gap                     <- D_i_s[2,"Y"] - D_i_s[1,"Y"]
				for(l in 1:(nrow(D_i_s)-1)){
                                        is_white_gap <- (D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap)
					if((s == "median" || s == "avg") & !is_white_gap )
						arrows(D_i_s[l,1],D_i_s[l,2],D_i_s[l+1,1],D_i_s[l+1,2],col=color_sel[["cluster_color"]][g],length=0,lwd=3)
                                        if((s == "lowquant" || s == "upquant") & !is_white_gap)
						arrows(D_i_s[l,1],D_i_s[l,2],D_i_s[l+1,1],D_i_s[l+1,2],col=color_sel[["cluster_color"]][g],length=0,lty="dashed",lwd=1)
                                        # else, as (is_white_gap == TRUE) then do not draw any arrow...
				}
				if(v == graphics$legend_on){
					# PLOT THE LEGEND ON ONE OF THE VIGNETTES			
					# DO THE LEGEND WITHIN [3.5;6.5] AND START FROM THE TOP, I.E. 6.5
					#y1 <- y0 <- 0.9*graphics[["ylim"]]$max-(g-1)*3/nrow(data[[1]])
					y1 <- y0 <- graphics[["legend_ycoord"]] + 5*g/nrow(data[[1]])
					x0 <- graphics[["legend_xcoord"]]
					x1 <- 1.2*graphics[["legend_xcoord"]]
					arrows(x0,y0,x1,y1,col=color_sel[["cluster_color"]][g],length=0,lwd=3)
                                        local_g <- color_sel[["cluster_order"]][g]
					text(x0+0.15,y1,labels=paste2(local_g," (",class_count[as.numeric(local_g)],")"),pos=4)
				}						
			}
			# DIMNAMES
			axis(2,at=D_i_s[,2],labels=dimnames(D_i_s)[[1]],las=2,tick=FALSE)
			axis(1,at=seq(from = graphics[["xlim"]]$min , to = graphics[["xlim"]]$max , by = ((graphics[["xlim"]]$max - graphics[["xlim"]]$min)/4)))
		}
	}
	# DEFINE COLOR GRADIENT
	gradient.colors <- graphics[["color_gradient"]] 
	# DETERMINES ROW AND COLUMN ORDERINGS OF THE AVG PATTERN
        data <- data[[1]][,canalysis_outcomes]
	rowD <- hclust(dist(data))
	colD <- hclust(dist(t(data)))
	# d <- data[rowD$order,colD$order]
        d <- data[rowD$order, rev(sort(colnames(data)))]
	# PRODUCE THE 3-PLOTS: IMAGE/HEATMAP, ROW AND COLUMN DENDROGRAMS
	image(1L:nrow(d), 1L:ncol(d), d, xlim = 0.5 + c(0, nrow(d)), ylim = 0.5 + c(0, ncol(d)), axes = FALSE, xlab = "", ylab = "",col=gradient.colors,main="Average pattern visualization")
	axis(2, 1L:ncol(d), labels = dimnames(d)[[2]], las = 2, line = -0.5, tick = 0, cex.axis = graphics[["cex"]])
	axis(1, 1L:nrow(d), labels = dimnames(d)[[1]], las = 1, line = -0.5, tick = 0, cex.axis = graphics[["cex"]])
	plot(rowD, axes = FALSE, yaxs = "i", main="Patterns arranged by similarity",ylab=NULL,cex=graphics[["cex"]])
	plot(colD, axes = FALSE, yaxs = "i", main="Outcomes arranged by similarity",ylab=NULL,cex=graphics[["cex"]])
}
# 2008-02-18: Exploratory Data Analysis for longitudinal data, profile plots
#
#
graphic_profile <- function(data,psoutput=NULL,plot_sample=20){
   if(!is.null(psoutput))
      postscript(file=psoutput,paper="a4")
   par(mfrow=c(3,3), 'las'=1,'mai'=c(0.3,0.3,0.3,0.5))
   id_set    <- sample(row.names(data))
   #
   for(o in colnames(data[,,1])){
      a <- as.data.frame(data[,o,])
      b <- reshape(a,direction="long",ids=row.names(a),times=1:dim(data)[[3]],timevar="Year",varying=list(names(a)),v.names=o)
      # INIT PLOT
      plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(min(a,na.rm=TRUE),max(a,na.rm=TRUE)),ylab=o)
      # RAW DATA
      for(i in 1:length(id_set)){
         if((length(id_set)-i) < plot_sample)
            line_col <- "blue"
         else
            line_col <- "black"
         lines(xy.coords(b[b$id == id_set[i],c("Year",o)]),type="l",col=line_col)
         }
      # STATS
      for(s in list(list(fun=function(x){return(mean(x,na.rm=TRUE))},col="red",lty="solid",lwd = 2),
                  list(fun=function(x){return(median(x,na.rm=TRUE))},col="red",lty="dashed",lwd = 1),
                  list(fun=function(x){return(quantile(x, na.rm=TRUE, probs = seq(0, 1, 0.05))["5%"])},col="red",lty="dashed",lwd = 1),
                  list(fun=function(x){return(quantile(x, na.rm=TRUE, probs = seq(0, 1, 0.05))["95%"])},col="red",lty="dashed",lwd = 1)
                  )){
         s[["pattern"]] <- cbind(Year=1:dim(data)[[3]],o=apply(data[,o,],2,s[["fun"]]))
         lines(xy.coords(s[["pattern"]]),type="l",col = s[["col"]],lty = s[["lty"]],lwd = s[["lwd"]])
      }
   }
   if(!is.null(psoutput))
      graphics.off()
}
# 2008-02-18: Exploratory Data Analysis for longitudinal data, profile plots
#
#
graphic_cluster_profiles <- function(data,canalysis_outcomes = parkinson_canalysis_outcomes,
        psoutput        =NULL,
        fun_stats       =list(
#                f_median   = list(fun = function(x){return(median(x,na.rm=TRUE))},lty="solid",lwd = 4)
                 f_mean   = list(fun = function(x){return(mean(x,na.rm=TRUE))},  lty="solid",lwd = 4)
                , f_25quant = list(fun = function(x){return(quantile(x,probs=0.25,na.rm=TRUE))},  lty="dashed",lwd = 1)
                , f_75quant = list(fun = function(x){return(quantile(x,probs=0.75,na.rm=TRUE))},  lty="dashed",lwd = 1)
                ))
                {
   #
   if(!is.null(psoutput))
      postscript(file=psoutput,paper="a4")
   #
   par(mfrow=c(2,9), 'las'=1,'mai'=c(0.3,0.1,0.3,0.3),oma=c(1,1,1,2))
   avg_pattern  <- compute_pattern(as.data.frame(data[,c(canalysis_outcomes,"class"),1]),fun = fun_stats[[1]][[1]])
   color_sel    <- get_coloring_scheme(avg_pattern,canalysis_outcomes = canalysis_outcomes)
   ymin         <- min(quantile(data,probs=0.05,na.rm=TRUE),na.rm=TRUE)
   ymax         <- max(quantile(data,probs=0.95,na.rm=TRUE),na.rm=TRUE)
   #
   id_set    <- sample(row.names(data))
   for(o in  canalysis_outcomes){
      # INIT PLOT
      plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(ymin,ymax),main=o,ann=TRUE,axes=FALSE)
      # STATS
      for(cl in 1:length(color_sel[["cluster_order"]])){
         df_subset <- data[data[,"class",1] == color_sel[["cluster_order"]][cl],,]
         df_subset <- df_subset[row.names(df_subset)[!is.na(row.names(df_subset))],,]
         # RAW DATA
         a <- as.data.frame(df_subset[,o,])
         b <- reshape(a,direction="long",ids=row.names(a),times=1:dim(data)[[3]],timevar="Year",varying=list(names(a)),v.names=o)
         # RAW PROFILES
#         for(i in sample(unique(b$id))[1:5])
#            lines(xy.coords(b[b$id == i,c("Year",o)]),type="l",col=color_sel[["cluster_color"]][cl],lwd=1,lty="dotted")
         for(s in fun_stats){
            s[["pattern"]] <- cbind(Year=1:dim(data)[[3]],o=apply(df_subset[,o,],2,s[["fun"]]))
            lines(xy.coords(s[["pattern"]]),type="l",col = color_sel[["cluster_color"]][cl],lty = s[["lty"]],lwd = s[["lwd"]])
         }
         # PLOT THE LEGEND ON ONE OF THE VIGNETTES			
         # DO THE LEGEND WITHIN [3.5;6.5] AND START FROM THE TOP, I.E. 6.5
      }
      # Y-AXIS PLOT (2)
      axis(2,at=pretty(c(ymin,ymax),n=3))
   }
   #
   # LEGEND
   #
   plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(ymin,ymax),main="legend",ann=TRUE,axes=FALSE)
   for(cl in 1:length(color_sel[["cluster_order"]])){
      yrange <- (ymax-ymin)/length(color_sel[["cluster_order"]])
      y1 <- y0 <- ymax-cl*0.5*yrange
      x0 <- 1
      x1 <- 2
      arrows(x0,y0,x1,y1,col=color_sel[["cluster_color"]][cl],length=0,lwd=3)
      text(1.15*x1,y1,labels=paste2(cl," (",nrow(data[data[,"class",1] == color_sel[["cluster_order"]][cl],,]),")"),pos=4)
   }
   #
   #
   #
   if(!is.null(psoutput))
      graphics.off()
}
# 2008-02-25:
#
# TAKE AS INPUT A DATA ARRAY
get_coloring_scheme <- function(data, canalysis_outcomes = parkinson_canalysis_outcomes){
   data <- data.matrix(data)
   # MAKE A DENDROGRAM FROM THESE AVG PATTERN FOR THE DENDROGRAM ORDERING
   dendro_patterns      <- hclust(dist(data))
   # SELECT DIVERGING COLORS FROM RColorBrewer
   cluster_colors       <- brewer.pal(nrow(data),"Set1")
   # RETURN A DATA STRUCTURE WITH CLUSTERS RIGHTLY ORDERED AND THEIR ASSOCIATED COLOR
   return(list(cluster_order=dendro_patterns$order,cluster_color=cluster_colors))
}
# 2007-07-03:
# SELECT ONLY ONE PAIR OF INDIVIDUALS IN EACH FAMILY. ADDITIONAL INDIVIDUALS
# IN LARGE FAMILIES WITH 3 OR 4 ASCERTAINED SIBLINGS ARE DISCARDED. FAMILIES
# WITH ONLY ONE INDIVIDUAL ARE DISCARDED. 
sibships_by_pairs	<- function(data){
	# REMOVE INDIVIDUALS IN SIBSHIPS >= 3
	individual_list <- dimnames(data[data[,"family"]  %in% names(table(data[,"family"])[table(data[,"family"])>=2]) ,])[[1]]
	# FOR EACH FAMILY (ROW), LIST THE ASCERTAINED SIBLINGS (COLUMNS: 1,2,3,4,5)
	sibling_list		<- reshape(data[individual_list,c("UNIEK","family","member")],
							timevar="member", idvar="family", direction="wide")
	siblingpair_list	<- c()
	# FOR EACH FAMILY, CONSIDER ONLY THE TWO FIRST SIBLINGS
	for(s in 1:nrow(sibling_list))
		siblingpair_list<- rbind(siblingpair_list,t(sibling_list[s,!is.na(sibling_list[s,])][2:3]))
	return(list(data=data[siblingpair_list,],model=NULL,type="sibpairs"))
}
# 2007-11-27: 
#
# L1, COLUMN NORMALIZATION, MAKE ALL OUTCOMES HAVE 1 AS MAXIMAL VALUE
transform_L1 <- function(data){
	data <- data/matrix(apply(data,2,sum,na.rm=TRUE),nrow(data),ncol(data),byrow=TRUE)
	# TODO: UPDATE TYPE FIELD	
	return(list(data=data,model=NULL,type="L1"))
}
# 2007-11-27: 
#
# L2, COLUMN NORMALIZATION, WEIGHT ALL OUTCOMES GIVEN THE SQUARE 'LENGTH'
transform_L2 <- function(data){
	data <- data/sqrt(matrix(apply(data^2,2,sum,na.rm=TRUE),nrow(data),ncol(data),byrow=TRUE))
	# TODO: UPDATE TYPE FIELD
	return(list(data=data,model=NULL,type="L2"))
}
# 2007-11-27: 
#
# SHIFT BACK TO ZERO ALL MEASUREMENTS
transform_MEDIAN <- function(data,model=NULL){
	return(transform(data=data,model=model,type="median"))
}
# 2007-11-27: 
#
# SHIFT BACK TO ZERO ALL MEASUREMENTS
transform_MIN <- function(data,model=NULL){
	return(transform(data=data,model=model,type="min"))
}
# 2007-11-27: 
#
# MAX
transform_ABSMAX <- function(data,model=NULL){
	return(transform(data=data,model=model,type="absmax"))
}
# 2007-11-27: 
#
# MAX
transform_AVG <- function(data,model=NULL){
	return(transform(data=data,model=model,type="avg"))
}
# 2007-11-27: 
#
# MAX
transform_MAX <- function(data,model=NULL){
	return(transform(data=data,model=model,type="max"))
}
# 2007-12-06:
#
#
transform <- function(data,model=NULL,type="min",effect="duration"){
	if(is.null(model)){
		run.vector <- dimnames(data)[[2]]
		model <- list()
		for(i in run.vector){
			if(type == "max"){
				model[[i]] <- max(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "sigma"){
				model[[i]] <- sd(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "absmax"){
				model[[i]] <- max(abs(data[,i],na.rm=TRUE),na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "min"){
				model[[i]] <- min(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			if(type == "avg"){
				model[[i]] <- mean(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			if(type == "median"){
				model[[i]] <- median(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			}
		return(list(data=data,type=type,model=model))
	}
	else{
		for(o in names(model)){
			if(type == "max" || type == "absmax" || type == "sigma")
				data[,o] <- data[,o]/model[[o]]
			if(type == "min" || type == "median" || type == "avg")
				data[,o] <- data[,o]-model[[o]]				
		}
		return(list(data=data,type=type,model=model))
	}
}
# 2008-02-18: 
#
# SIGMA
transform_SIGMA <- function(data,model=NULL){
	return(transform(data=data,model=model,type="sigma"))
}
# 2008-02-26:
#
#
transform_longitudinal  <- function(data, fun, canalysis_outcomes = parkinson_canalysis_outcomes,effect="duration"){
   data_out             <- data
   m                    <- fun(as.data.frame(data[,c(canalysis_outcomes,effect),1]))
   m_effect             <- 
   data_out[row.names(m[["data"]]),colnames(m[["data"]]),1]     <- data.matrix(m[["data"]])
   for(i in 2:dim(data)[[3]]){
      # WE WANT TO REMOVE THE EFFECT OF DISEASE DURATION (AKA THE TIME) FOR BASELINE,
      # BUT WE PRESERVE THE TIME EFFECT IN THE LONGITUDINAL ANALYSIS. 
      data_tmp          <- as.data.frame(cbind(data[,canalysis_outcomes,i],data[,effect,1]))
      colnames(data_tmp)[ncol(data_tmp)] <- effect 
      m_tmp          <- fun(data_tmp,model = m[["model"]])
      data_out[row.names(m_tmp[["data"]]),canalysis_outcomes,i]  <- data.matrix(m_tmp[["data"]])[,canalysis_outcomes]
      }
   return(data_out)
}
# 2008-03-02
#
# PROBANT <- 1 ; SIBLING <- 2
transform_sibordering <- function(data){
   for(f in unique(data$family)){
      fmembers <- sort(data[data$family == f,"member"])
      if(!is.na(fmembers[1]) && !is.na(fmembers[2])){
         data[data$family == f & data$member == fmembers[1],"member"] <- 1
         data[data$family == f & data$member == fmembers[2],"member"] <- 2 
      }
   }
   return(data[data$member <= 2,])
}
# 2008-03-20:
#
#
transform_adjust <- function(
        data,
        model=NULL,
        type="lm",
        f_matrix=formula_matrix,
        transform=list(
                 sqrt=  list(ls="sqrt(", rs=")")
                ,log=   list(ls="log(",  rs=")")
                ,none=  list(ls="",      rs="")
                ,square=list(ls="(",     rs=")^2")
                ,exp=   list(ls="exp(",  rs=")")
                ),
        optimize = FALSE
        )
{
   if(is.null(model)){
      # data
      data_out     <- as.data.frame(matrix(NA,nrow(data),0,dimnames=list(row.names(data),list())))
      model_out    <- list()
      # for each outcome (column of the data matrix), loop
      for(i in 1:nrow(f_matrix)){
         fact_name <- colnames(f_matrix)[i]
         # left hand side of the formula (above 0)
         lhs       <- colnames(f_matrix)[f_matrix[i,]>0] 
         # right hand side of the formula (below 0)
         rhs       <- colnames(f_matrix)[f_matrix[i,]<0]
         if(length(lhs) == 1 && length(rhs) > 0){
            list_f    <- list(paste2(paste(lhs,collapse="+"),"~"))
            # for each left hand side element
            for(f_idx in 1:length(list_f)){
               elmt <- list_f[[f_idx]]
               # we append the right hand side transform in only one way for both
               for(t_id in 1:length(transform)){
                  # in case this is the first transform, we just update list_f,
                  # otherwise we add new elements to list_f referred to as f_idx_l
                  f_idx_l <- f_idx+t_id-1
                  local_f <- ""
                  for(rhs_l in rhs){
                     local_f <- paste2(local_f,transform[[t_id]]$ls,rhs_l,transform[[t_id]]$rs)
                     if(length(rhs) > 1 && rhs_l != rhs[[length(rhs)]])
                        local_f <- paste2(local_f,"+")
                  }
                  list_f[[f_idx_l]] <- paste2(elmt,local_f)
               }
            }
            model          <- list()
            adj_rsquare    <- matrix(0,length(list_f))
            for(f_id in 1:length(list_f)){
               model[[f_id]]          <- lm(as.formula(list_f[[f_id]]),data=data)
               adj_rsquare[f_id]      <- summary(model[[f_id]])[["adj.r.squared"]]
            }
            # if there are several rsquare that are equal, take the 1st one...
            best_id                            <- which(adj_rsquare == max(adj_rsquare))[1]
            model_out[[fact_name]] <- model[[best_id]]
            model_out[[fact_name]][["rhs"]]      <- rhs
            model_out[[fact_name]][["formula"]] <- list_f[[best_id]]
            # take the corresponding residual
            r                              <- residuals(model_out[[fact_name]])
            r_names                        <- names(r)
            }# else, do not adjust data and preserve original scores
         else
            r <- NULL
         data_out                       <- cbind(data_out,matrix(NA,nrow(data_out),1))
         if(!is.null(r))
            data_out[r_names,ncol(data_out)] <- r
         else
            data_out[,ncol(data_out)]   <- data[row.names(data_out),fact_name]
         colnames(data_out)[i]          <- colnames(f_matrix)[i]
      }
      return(list(data=data_out,model=model_out,type="lm"))
   }
   else{
      for(o in names(model)){
         data[,o] <- data[,o]-predict(model[[o]],data[,c(o,unlist(model[[o]][["rhs"]]))])
         }
      return(list(data=data,type=type,model=model))
   }
}
# 2008-03-29:
#
#
transform_addnoise <- function(
        data,
        model=NULL,
        canalysis_outcomes = NULL,
        type="addnoise",
        relative_noise_amount = 0.01,
        rand_seed = 6013
        )
{
   if(is.null(model)){
      for(o in which(canalysis_outcomes %in% colnames(data))){
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
   return(list(data=data,type=type,model=model))
}
# 2007-07-03:
# SELECT ONLY ONE PAIR OF INDIVIDUALS IN EACH FAMILY. ADDITIONAL INDIVIDUALS
# IN LARGE FAMILIES WITH 3 OR 4 ASCERTAINED SIBLINGS ARE DISCARDED. FAMILIES
# WITH ONLY ONE INDIVIDUAL ARE DISCARDED. 
sibships_by_pairs	<- function(data){
	# REMOVE INDIVIDUALS IN SIBSHIPS >= 3
	individual_list <- dimnames(data[data[,"family"]  %in% names(table(data[,"family"])[table(data[,"family"])>=2]) ,])[[1]]
	# FOR EACH FAMILY (ROW), LIST THE ASCERTAINED SIBLINGS (COLUMNS: 1,2,3,4,5)
	sibling_list		<- reshape(data[individual_list,c("UNIEK","family","member")],
							timevar="member", idvar="family", direction="wide")
	siblingpair_list	<- c()
	# FOR EACH FAMILY, CONSIDER ONLY THE TWO FIRST SIBLINGS
	for(s in 1:nrow(sibling_list))
		siblingpair_list<- rbind(siblingpair_list,t(sibling_list[s,!is.na(sibling_list[s,])][2:3]))
	return(list(data=data[siblingpair_list,],model=NULL,type="sibpairs"))
}
# 2007-11-27: 
#
# L1, COLUMN NORMALIZATION, MAKE ALL OUTCOMES HAVE 1 AS MAXIMAL VALUE
transform_L1 <- function(data){
	data <- data/matrix(apply(data,2,sum,na.rm=TRUE),nrow(data),ncol(data),byrow=TRUE)
	# TODO: UPDATE TYPE FIELD	
	return(list(data=data,model=NULL,type="L1"))
}
# 2007-11-27: 
#
# L2, COLUMN NORMALIZATION, WEIGHT ALL OUTCOMES GIVEN THE SQUARE 'LENGTH'
transform_L2 <- function(data){
	data <- data/sqrt(matrix(apply(data^2,2,sum,na.rm=TRUE),nrow(data),ncol(data),byrow=TRUE))
	# TODO: UPDATE TYPE FIELD
	return(list(data=data,model=NULL,type="L2"))
}
# 2007-11-27: 
#
# SHIFT BACK TO ZERO ALL MEASUREMENTS
transform_MEDIAN <- function(data,model=NULL){
	return(transform(data=data,model=model,type="median"))
}
# 2007-11-27: 
#
# SHIFT BACK TO ZERO ALL MEASUREMENTS
transform_MIN <- function(data,model=NULL){
	return(transform(data=data,model=model,type="min"))
}
# 2007-11-27: 
#
# MAX
transform_ABSMAX <- function(data,model=NULL){
	return(transform(data=data,model=model,type="absmax"))
}
# 2007-11-27: 
#
# MAX
transform_AVG <- function(data,model=NULL){
	return(transform(data=data,model=model,type="avg"))
}
# 2007-11-27: 
#
# MAX
transform_MAX <- function(data,model=NULL){
	return(transform(data=data,model=model,type="max"))
}
# 2007-12-06:
#
#
transform <- function(data,model=NULL,type="min",effect="duration"){
	if(is.null(model)){
		run.vector <- dimnames(data)[[2]]
		model <- list()
		for(i in run.vector){
			if(type == "max"){
				model[[i]] <- max(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "sigma"){
				model[[i]] <- sd(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "absmax"){
				model[[i]] <- max(abs(data[,i],na.rm=TRUE),na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "min"){
				model[[i]] <- min(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			if(type == "avg"){
				model[[i]] <- mean(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			if(type == "median"){
				model[[i]] <- median(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			}
		return(list(data=data,type=type,model=model))
	}
	else{
		for(o in names(model)){
			if(type == "max" || type == "absmax" || type == "sigma")
				data[,o] <- data[,o]/model[[o]]
			if(type == "min" || type == "median" || type == "avg")
				data[,o] <- data[,o]-model[[o]]				
		}
		return(list(data=data,type=type,model=model))
	}
}
# 2008-02-18: 
#
# SIGMA
transform_SIGMA <- function(data,model=NULL){
	return(transform(data=data,model=model,type="sigma"))
}
# 2008-02-26:
#
#
transform_longitudinal  <- function(data, fun, canalysis_outcomes = parkinson_canalysis_outcomes,effect="duration"){
   data_out             <- data
   m                    <- fun(as.data.frame(data[,c(canalysis_outcomes,effect),1]))
   m_effect             <- 
   data_out[row.names(m[["data"]]),colnames(m[["data"]]),1]     <- data.matrix(m[["data"]])
   for(i in 2:dim(data)[[3]]){
      # WE WANT TO REMOVE THE EFFECT OF DISEASE DURATION (AKA THE TIME) FOR BASELINE,
      # BUT WE PRESERVE THE TIME EFFECT IN THE LONGITUDINAL ANALYSIS. 
      data_tmp          <- as.data.frame(cbind(data[,canalysis_outcomes,i],data[,effect,1]))
      colnames(data_tmp)[ncol(data_tmp)] <- effect 
      m_tmp          <- fun(data_tmp,model = m[["model"]])
      data_out[row.names(m_tmp[["data"]]),canalysis_outcomes,i]  <- data.matrix(m_tmp[["data"]])[,canalysis_outcomes]
      }
   return(data_out)
}
# 2008-03-02
#
# PROBANT <- 1 ; SIBLING <- 2
transform_sibordering <- function(data){
   for(f in unique(data$family)){
      fmembers <- sort(data[data$family == f,"member"])
      if(!is.na(fmembers[1]) && !is.na(fmembers[2])){
         data[data$family == f & data$member == fmembers[1],"member"] <- 1
         data[data$family == f & data$member == fmembers[2],"member"] <- 2 
      }
   }
   return(data[data$member <= 2,])
}
# 2008-03-20:
#
#
transform_adjust <- function(
        data,
        model=NULL,
        type="lm",
        f_matrix=formula_matrix,
        transform=list(
                 sqrt=  list(ls="sqrt(", rs=")")
                ,log=   list(ls="log(",  rs=")")
                ,none=  list(ls="",      rs="")
                ,square=list(ls="(",     rs=")^2")
                ,exp=   list(ls="exp(",  rs=")")
                ),
        optimize = FALSE
        )
{
   if(is.null(model)){
      # data
      data_out     <- as.data.frame(matrix(NA,nrow(data),0,dimnames=list(row.names(data),list())))
      model_out    <- list()
      # for each outcome (column of the data matrix), loop
      for(i in 1:nrow(f_matrix)){
         fact_name <- colnames(f_matrix)[i]
         # left hand side of the formula (above 0)
         lhs       <- colnames(f_matrix)[f_matrix[i,]>0] 
         # right hand side of the formula (below 0)
         rhs       <- colnames(f_matrix)[f_matrix[i,]<0]
         if(length(lhs) == 1 && length(rhs) > 0){
            list_f    <- list(paste2(paste(lhs,collapse="+"),"~"))
            # for each left hand side element
            for(f_idx in 1:length(list_f)){
               elmt <- list_f[[f_idx]]
               # we append the right hand side transform in only one way for both
               for(t_id in 1:length(transform)){
                  # in case this is the first transform, we just update list_f,
                  # otherwise we add new elements to list_f referred to as f_idx_l
                  f_idx_l <- f_idx+t_id-1
                  local_f <- ""
                  for(rhs_l in rhs){
                     local_f <- paste2(local_f,transform[[t_id]]$ls,rhs_l,transform[[t_id]]$rs)
                     if(length(rhs) > 1 && rhs_l != rhs[[length(rhs)]])
                        local_f <- paste2(local_f,"+")
                  }
                  list_f[[f_idx_l]] <- paste2(elmt,local_f)
               }
            }
            model          <- list()
            adj_rsquare    <- matrix(0,length(list_f))
            for(f_id in 1:length(list_f)){
               model[[f_id]]          <- lm(as.formula(list_f[[f_id]]),data=data)
               adj_rsquare[f_id]      <- summary(model[[f_id]])[["adj.r.squared"]]
            }
            # if there are several rsquare that are equal, take the 1st one...
            best_id                            <- which(adj_rsquare == max(adj_rsquare))[1]
            model_out[[fact_name]] <- model[[best_id]]
            model_out[[fact_name]][["rhs"]]      <- rhs
            model_out[[fact_name]][["formula"]] <- list_f[[best_id]]
            # take the corresponding residual
            r                              <- residuals(model_out[[fact_name]])
            r_names                        <- names(r)
            }# else, do not adjust data and preserve original scores
         else
            r <- NULL
         data_out                       <- cbind(data_out,matrix(NA,nrow(data_out),1))
         if(!is.null(r))
            data_out[r_names,ncol(data_out)] <- r
         else
            data_out[,ncol(data_out)]   <- data[row.names(data_out),fact_name]
         colnames(data_out)[i]          <- colnames(f_matrix)[i]
      }
      return(list(data=data_out,model=model_out,type="lm"))
   }
   else{
      for(o in names(model)){
         data[,o] <- data[,o]-predict(model[[o]],data[,c(o,unlist(model[[o]][["rhs"]]))])
         }
      return(list(data=data,type=type,model=model))
   }
}
# 2008-03-29:
#
#
transform_addnoise <- function(
        data,
        model=NULL,
        canalysis_outcomes = NULL,
        type="addnoise",
        relative_noise_amount = 0.01,
        rand_seed = 6013
        )
{
   if(is.null(model)){
      for(o in which(canalysis_outcomes %in% colnames(data))){
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
   return(list(data=data,type=type,model=model))
}
# 2006-12-11:
#
#
mclust2   <- function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = NULL, warn = FALSE, ...) 
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
# 2006-12-11:
#
#
summary_mclust2_BIC <- function (object, data, G = NULL, modelNames = NULL, ...) 
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
# 2006-12-11:
#
#
mclust2_bicaic  <- function (modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE,...)
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
# 2006-12-11:
# 
#
mclust2_BIC        <- function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), 
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
#
# ADDITIONAL DATA DEFINITION
#
f_sqrt  <- list(ls="sqrt(", rs=")")
f_log   <- list(ls="log(",  rs=")")
f_none  <- list(ls="",      rs="")
f_square<- list(ls="(",     rs=")^2")
f_exp   <- list(ls="exp(",  rs=")")
#
# PARKINSON CONFIG DATA
#
#
parkinson_canalysis_outcomes <- sort(c("spestremor", "spesbrady", "spesrigiditeit", "spesaxial",
   "cogmemory", "cogattention", "cogexecutief", "cogvisuo",
   "autgi", "auturin", "autcardio",
   "motfluct", "dysk", "psychosis",
   "sleep", "eds", "beck"))
#                        "duration", # NOT INVOLVED IN THE CLUSTER ANALYSIS
parkinson_stat_groups <- list(Motor=c("spestremor","spesbrady","spesrigiditeit","spesaxial"),
   Cognition=c("cogmemory","cogattention","cogexecutief","cogvisuo"),
   Autonomic=c("autgi","auturin","autcardio"),
   Motfluct=c("motfluct"),
   Dyskinesia=c("dysk"),
   Psychosis=c("psychosis"),
   Sleep=c("sleep"),
   Eds=c("eds"),
   Beck=c("beck"),
   Duration=c("duration"))
#
parkinson_visu_matrix <- matrix(c(
"cogmemory",0,          "9,5","Cognitive and Motricity" #1
,"cogattention",0,       "6,5","Cognitive and Motricity"
,"cogexecutief",0,       "8,5","Cognitive and Motricity"
,"cogvisuo",0,           "7,5","Cognitive and Motricity"
,"spestremor",0,         "4,5","Cognitive and Motricity" #5
,"spesbrady",0,          "2,5","Cognitive and Motricity"
,"spesrigiditeit",0,     "3,5","Cognitive and Motricity"
,"spesaxial",0,          "1,5","Cognitive and Motricity"
,"autgi",0,              "8,5","Autonomic and Sumscores"
,"auturin",0,            "9,5","Autonomic and Sumscores" #10
,"autcardio",0,          "7,5","Autonomic and Sumscores"
,"motfluct",0,           "4,5","Autonomic and Sumscores"
,"dysk",0,               "3,5","Autonomic and Sumscores"
,"psychosis",0,          "8,5","Other"
,"sleep",0,              "5,5","Other" #15
,"eds",0,                "6,5","Other"
,"beck",0,               "7,5","Other"
,"duration",0,           "9,5","Other"),18,4,dimnames=list(list(),list("Pheno","X","Y","Z")))
#
# GARP CONFIG DATA
# 
garp_canalysis_outcomes <- c("HP_L","HP_R","KN_L","KN_R","DIP2_L","DIP3_L","DIP4_L","DIP5_L"
   ,"PIP2_L","PIP3_L","PIP4_L","PIP5_L","IP_L","CMC_L","DIP2_R","DIP3_R"
   ,"DIP4_R","DIP5_R","PIP2_R","PIP3_R","PIP4_R","PIP5_R","IP_R","CMC_R"
   ,"CSF12","CSF23","CSF34","CSF45","CSF56","CSF67","LSF12","LSF23"
   ,"LSF34","LSF45","LSF56","LSD12","LSD23","LSD34","LSD45","LSD56"
   ,"CSD23","CSD34","CSD45","CSD56","CSD67")
garp_stat_groups <- list(
   demographic = c("age","bmi")
   ,hips = c("HP_L","HP_R")
   ,knees = c("KN_L","KN_R")
   ,cmc = c("CMC_L","CMC_R")
   ,dips = c("IP_L","DIP2_L","DIP3_L","DIP4_L","DIP5_L","IP_R","DIP2_R","DIP3_R","DIP4_R","DIP5_R")
   ,pips = c("PIP2_L","PIP3_L","PIP4_L","PIP5_L","PIP2_R","PIP3_R","PIP4_R","PIP5_R")
   ,facets = c("CSF12","CSF23","CSF34","CSF45","CSF56","CSF67","LSF12","LSF23","LSF34","LSF45","LSF56")
   ,discus = c("LSD12","LSD23","LSD34","LSD45","LSD56","CSD23","CSD34","CSD45","CSD56","CSD67"))
garp_visu_matrix <- matrix(c(
"fac1_4",0,10,"Factor Analysis and Demo data"#1
,"fac2_4",0,9,"Factor Analysis and Demo data"
,"fac3_4",0,8,"Factor Analysis and Demo data"
,"fac1_1",0,6,"Factor Analysis and Demo data"
,"fac2_1",0,5,"Factor Analysis and Demo data"#5
,"fac3_1",0,4,"Factor Analysis and Demo data"
,"age",0,2,"Factor Analysis and Demo data"
,"bmi",0,1,"Factor Analysis and Demo data"
,"HP_L",0,10,"Hips, Knees, CMC"
,"HP_R",0,9,"Hips, Knees, CMC"#10
,"KN_L",0,6,"Hips, Knees, CMC"
,"KN_R",0,5,"Hips, Knees, CMC"
,"CMC_L",0,2,"Hips, Knees, CMC"
,"CMC_R",0,1,"Hips, Knees, CMC"
,"DIP2_L",0,8,"Hands DIP_IP"#15
,"DIP3_L",0,9,"Hands DIP_IP"
,"DIP4_L",0,10,"Hands DIP_IP"
,"DIP5_L",0,11,"Hands DIP_IP"
,"PIP2_L",0,7,"Hands PIP"
,"PIP3_L",0,8,"Hands PIP"#20
,"PIP4_L",0,9,"Hands PIP"
,"PIP5_L",0,10,"Hands PIP"
,"IP_L",0,7,"Hands DIP_IP"
,"DIP2_R",0,4,"Hands DIP_IP"
,"DIP3_R",0,3,"Hands DIP_IP"#25
,"DIP4_R",0,2,"Hands DIP_IP"
,"DIP5_R",0,1,"Hands DIP_IP"
,"PIP2_R",0,4,"Hands PIP"
,"PIP3_R",0,3,"Hands PIP"
,"PIP4_R",0,2,"Hands PIP"#30
,"PIP5_R",0,1,"Hands PIP"
,"IP_R",0,5,"Hands DIP_IP"
,"CSF12",0,11,"Spine Facets"
,"CSF23",0,10,"Spine Facets"
,"CSF34",0,9,"Spine Facets"#35
,"CSF45",0,8,"Spine Facets"
,"CSF56",0,7,"Spine Facets"
,"CSF67",0,6,"Spine Facets"
,"LSF12",0,4,"Spine Facets"
,"LSF23",0,3,"Spine Facets"#40
,"LSF34",0,2,"Spine Facets"
,"LSF45",0,1,"Spine Facets"
,"LSF56",0,0,"Spine Facets"
,"LSD12",0,4,"Spine Discus"
,"LSD23",0,3,"Spine Discus"#45
,"LSD34",0,2,"Spine Discus"
,"LSD45",0,1,"Spine Discus"
,"LSD56",0,0,"Spine Discus"
,"CSD23",0,10,"Spine Discus"
,"CSD34",0,9,"Spine Discus"#50
,"CSD45",0,8,"Spine Discus"
,"CSD56",0,7,"Spine Discus"
,"CSD67",0,6,"Spine Discus"),53,4,dimnames=list(list(),list("Pheno","X","Y","Z")))
