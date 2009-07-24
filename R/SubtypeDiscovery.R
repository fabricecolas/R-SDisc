`analysis` <-
function(x, device="PS", img=FALSE, html=TRUE){
   plot(attr(x,"cdata"), device = device)
   x <- fun_cmodel(x)
   file_save <- paste2(attr(x,"prefix"),"_IMAGE.RData")
   cat("\nSave 'cresult' into ",file_save)
   save(list='x', file=file_save)
   cat("\nWrite best models\n")
   write.cresult(x)
   cat("\nPlot cresult -> PS")
   plot(x, device=device)
   cat("\nGenerate HTML report -> .html")
   print(x, html=html, img=img)

   return(x)
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



`stability_analysis` <-
function(x,q,rseed=6013,nnoise=10,nrep=10,ps=TRUE, ncolors=9,minmax=c(0.5,1)){
   settings <- strsplit(q,',')[[1]]
   cdata <- get_cdata(attr(x, "cdata"))
   r <- list()
   v <- matrix(NA,nnoise,nrep*(nrep-1)/2,dimnames=list(sprintf("%4f",(1/2)^(1:nnoise)),list()))
   # NOISE LEVELS
   for(sd_l in (1/2)^(1:nnoise)){ 
      l <- sprintf("%4f",sd_l)
      r[[l]] <- list()
      for(ridx in 1:nrep){
         s <- rseed+ridx
         r[[l]][[s]] <- list()
         set.seed(s)
         Y <- cdata+matrix(rnorm(nrow(cdata)*ncol(cdata),sd=sd_l,mean=0),nrow(cdata),ncol(cdata))
         r[[l]][[s]] <- fun_mbc_em(Y,modelName=settings[1],G=settings[2],rseed=as.numeric(settings[3]))
         r[[l]][[s]][["Y"]] <- Y
      }
      tmp <- c()
      for(ridx1 in 1:nrep){
         s1 <- rseed+ridx1
         for(ridx2 in 1:nrep){
            s2 <- rseed+ridx2
            if(ridx2>ridx1){
               m <- as.table(ftable(r[[l]][[s1]]$labelling,r[[l]][[s2]]$labelling))
               tmp <- c(tmp,sqrt(as.numeric(summary(m)["statistic"])/(as.numeric(summary(m)["n.cases"])*(min(dim(m))-1))))
            }
         }
      }
      v[l,] <- sort(tmp)
   }
   breaks <- seq(minmax[1],minmax[2],(minmax[2]-minmax[1])/ncolors)
   if(ps){
      postscript(paste2(attr(x,'prefix'),"_stability_analysis_",gsub(',','_',q),".ps"))
      image(v,col=brewer.pal(ncolors,"Greys"),breaks=breaks,axes=FALSE,ylab='Association level V, quantiles',cex.lab=1.7)
      contour(v,add=TRUE,labcex=1.3)
      axis(1,at=seq(0,1,0.11),labels=sprintf("%.1f%%",100*as.numeric(row.names(v))),las=2,cex.axis=1.7)
      axis(2,at=seq(0,1,0.2),labels=sprintf("%d%%",100*seq(0,1,0.2)),cex.axis=1.7)
      title(paste(settings,sep="",collapse=" "))
      graphics.off()
   }
   return(v)
}

`compare_transform` <-
function(x,y,probs=seq(0,1,0.1)){
   prefix               <- attr(x,"prefix")
   m <- data.frame(attr(x,'cfun_params'),V=NA)
   row.names(m) <- apply(as.matrix(m[,1:3]),1,paste,sep='',collapse=',') 
   for(i in names(x)){
      tmp <- table(data.frame(x=attr(x[[i]],'model')$labelling, y=attr(y[[i]],'model')$labelling))
      tmp_summary <- summary(tmp)
      n <- as.numeric(tmp_summary["n.cases"])
      X2 <- as.numeric(tmp_summary["statistic"])
      k <- min(dim(tmp))
      m[i,'V'] <- sqrt(X2/(n*(k-1)))
   }
   m <- na.omit(m)
   m_x <- xtabs(formula=V~modelName+G,aggregate(data.frame(V=m[,'V']),by=list(modelName=m[,'modelName'],G=m[,'G']),mean,na.rm=TRUE))
   return(list(summary=m_x, m=m))
}

`compare_cresult` <-
function(x,query=NULL){
   prefix               <- attr(x,"prefix")
   nbr_top_models       <- attr(x,"nbr_top_models")
   cross_comparison_results <- list()
   best_models          <- get_best_models(x)
   model_pairs          <- t(combn(best_models,m=2))
   if(nrow(model_pairs)>=1){
      cnames            <- gsub(",","_",paste("(",apply(model_pairs,1,paste,collapse="),("),")",sep=""))
      for(mp_idx in 1:nrow(model_pairs)){
         a              <- x[[which(model_pairs[mp_idx,1] == names(x))]]
         b              <- x[[which(model_pairs[mp_idx,2] == names(x))]]
         if( !is.na(attr(a,"model")$loglik) && !is.na(attr(b,"model")$loglik))
            cross_comparison_results[[cnames[mp_idx]]] <- cross_compare_models(a,b)
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
   m <- as.table(ftable(attr(x,"model")$labelling,attr(y,"model")$labelling))
   m_num <- m
   # USE JOINT PROBABILITIES TO COMPUTE THE CONTINGENCY TABLE (%) 
   m_p <- fullproba_ftable(attr(x,"model")$z,attr(y,"model")$z)
   map_rand_err_test <- t.test(as.numeric(m/sum(m)-m_p))
   # CHI2_STATS: TEST THE ASSOCIATION BETWEEN THE TWO CATEGORICAL VARIABLES
   # H0: NO ASSOCIATION BETWEEN THE TWO VARIABLES
   # H1: THERE IS ASSOCIATION BETWEEN THE TWO VARIABLES
   ind_test <- chisq.test(m,simulate.p.value=TRUE)
   # CRAMER'S V 
   x_summary <- summary(m)
   n <- as.numeric(x_summary["n.cases"]) 
   X2 <- as.numeric(x_summary["statistic"]) 
   k <- min(dim(m)) 
   V <- sqrt(X2/(n*(k-1))) 
   # PEARSON'S CORRELATION: REPORT CORRELATION BETWEEN THE TWO CATEGORICAL
   # VARIABLES TEST (P) WHETHER THE CORRELATION IS EQUAL TO ZERO AND REPORTS
   # THE 95% CONFIDENCE INTERVAL
   cor_test     <- cor.test(attr(x,"model")$labelling,attr(y,"model")$labelling,use="pairwise.complete.obs")
   # BIND STATISTICS SUCH AS LOG-ODD RATIO, LAMBDA-SIBS INTO A FINAL TABLE 1RST
   # DIRECTION
   m_html <- apply(m,1:2,html_format_numeric,fmt="%d")
   m[m_num == 0] <- m_html[m_num == 0] <- ""
   for(s in names(attr(x,"stats"))){
      m_html <- cbind(m_html,data.matrix(attr(x,"stats")[[s]][["out"]]))
      m	<- cbind(m,data.matrix(attr(x,"stats")[[s]][[2]]))
      }
   mat_stats <- mat_stats_html <- matrix(0,0,attr(y,"model")$G)
   for(s in names(attr(y,"stats"))){
      mat_stats_html <- rbind(mat_stats_html, t(attr(y,"stats")[[s]][["out"]]))
      mat_stats <- rbind(mat_stats, t(attr(y,"stats")[[s]][[2]]))
      }
   mat_stats <- cbind(mat_stats, matrix(0,nrow(mat_stats),ncol(m)-attr(y,"model")$G))
   mat_stats_html <- cbind(mat_stats_html, matrix(0,nrow(mat_stats),ncol(m)-attr(y,"model")$G))
   dimnames(mat_stats)[[2]] <- dimnames(mat_stats_html)[[2]] <-  dimnames(m)[[2]]
   m <- rbind(m,mat_stats)
   m_html <- rbind(m_html,mat_stats_html)
   # FINALLY, ORDER THE ROWS AND THE COLUMNS BY (DI)SIMILARIY, HIERARCHICAL
   # CLUSTERING
   ordering <- list(row=hclust(dist(m_num[1:attr(x,"model")$G,1:attr(y,"model")$G]))$order,
                     column=hclust(dist(t(m_num[1:attr(x,"model")$G,1:attr(y,"model")$G])))$order)
   m <- m[c(ordering$row,(attr(x,"model")$G+1):nrow(m)),
      c(ordering$col,(attr(y,"model")$G+1):ncol(m))]
   m_html <- m_html[c(ordering$row,(attr(x,"model")$G+1):nrow(m)),
      c(ordering$col,(attr(y,"model")$G+1):ncol(m))]
   # BIND CHI2 STATS
   m <- rbind(m,"")
   m_html <- rbind(m_html,"")
   row.names(m)[nrow(m)] <- row.names(m_html)[nrow(m_html)] <- 'Chi2'
   m[nrow(m),1:2] <- m_html[nrow(m_html),1:2] <- c(sprintf('%.1e',ind_test$p.value),sprintf('%.1f',ind_test$statistic))
   # BIND CRAMER'S V ASSOCIATION STAT
   m <- rbind(m,"")
   m_html <-  rbind(m_html,"")
   row.names(m)[nrow(m)] <- row.names(m_html)[nrow(m_html)] <- "Cramer's V"
   m_html[nrow(m_html),1] <- html_format_numeric(V,fmt='%.2f') 
   m[nrow(m),1] <- V
   # T-TEST TO EVALUATE WHETHER THE DIFFERENCES BETWEEN THE MAP CONTINGENCY
   # TABLE AND THE JOINT PROBABILITY DISTRIBUTION OF THE TWO FACTORS IS A
   # RANDOM ERROR OF MEAN = 0
   m <- rbind(m,"")
   m_html <-  rbind(m_html,"")
   row.names(m)[nrow(m)] <- row.names(m_html)[nrow(m_html)] <- "t-test (mapping)"
   m[nrow(m),1] <- m_html[nrow(m_html),1] <- sprintf("%.1e",map_rand_err_test$p.value)
   # TO QUICKEN THE LATER OPENING OF THE MANY CSV FILES,  WE REDUCE THE NUMBER
   # OF ZEROS TO CLEAR MANUALLY. IN WRITE.CSV2 'NA' ARE SUBSTITUTED BY EMPTY
   # STRINGS
   m[m == "0"] <- m_html[m_html == "0"] <- ""
   m <- as.data.frame(m)
   m_html <- as.data.frame(m_html)
   return(list(html=m_html,full=m))
}

`fun_mbc_em` <-
function(data, G=3, modelName="EII", rseed=6013) {
   G <- as.numeric(G)
   modelName <- as.character(modelName)
   # PREPARE A RANDOM-UNIFORM MATRIX OF CLUSTER MEMBERSHIPS WITH G COLS AND N
   # ROWS (MIN IS 0 AND MAX 1)
   if(G >= 2){
      set.seed(rseed)
      cmembership_matrix <- matrix(runif(G*nrow(data)), nrow(data), G)
      cmembership_matrix <- cmembership_matrix/apply(cmembership_matrix, 1, sum)
   }
   else
      cmembership_matrix <- matrix(1,nrow(data),G)
   # ESTIMATE AN INITIAL MODEL WHICH WE PROVIDE AS STARTING POINT TO EM
   mstep_est <- mstep(data=data, modelName=modelName, z=cmembership_matrix,
                  warn=FALSE)
   # AND THEN, STARTS EM WITH THE E-STEP
   model <- em(data=data, modelName=mstep_est$modelName,
               parameters=mstep_est$parameters, warn=FALSE)
   model[["labelling"]] <- map(model[["z"]],warn=FALSE)
   # CALCULATE THE BIC SCORE
   model[["BIC"]] <- bic(modelName=modelName, loglik=model$loglik, n=nrow(data)
                           , d=ncol(data), G=G)
   return(model)
}

`stats_generalization` <-
function(data,class
   , fun_classifier=list(model=model_naive_bayes, predict=predict_naive_bayes)
   , K=7, fun_eval=stratified_traintest_split, title = NULL) {
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
function(data, grBreak=6){
   # DEFINE THE DEFAULT VALUES OF OUR CONF MATRIX
   conf_col_names <- list("group","in_canalysis","fun_transform","visu_groups","visu_ycoord","heatmap_ycoord")
   default_values <- c("factor_1",TRUE,"transform_AVG transform_SIGMA","var_group_1","NA","NA")
   m <- matrix(default_values,ncol(data),length(default_values),byrow=TRUE,
                  dimnames=list(colnames(data),conf_col_names))
   m[,'heatmap_ycoord'] <- 1:nrow(m)
   mSeq <- 1+as.integer(nrow(m)/grBreak)
   mGroup <- matrix(c('group',NA,NA),nrow(m),3,byrow=TRUE)
   for(i in 0:(grBreak-1)){
      mGroup[(i*mSeq+1):min((i+1)*mSeq,nrow(m)),2] <- i+1
      mGroup[(i*mSeq+1):min((i+1)*mSeq,nrow(m)),3] <- 1:(min((i+1)*mSeq,nrow(m))-(i*mSeq))
   }
   m[,'group'] <- m[,'visu_groups'] <- apply(mGroup[,1:2],1,paste,collapse='_')
   m[,'visu_ycoord'] <- mGroup[,3]
   return(m)
}

`get_cdata_prcomp` <- function(x, cumproportion=0.95, scale=TRUE, center=FALSE){
   # RETRIEVE PREVIOUS PARAMETERS
   s1 <- as.matrix(attr(x,"settings"))
   prefix <- gsub('[0-9]{4}-[0-9]{2}-[0-9]{2}_','',paste2(attr(x,"prefix"),"_ON_PRCOMP_"))
   # PRINCIPAL COMPONENTS, 95%
   pcomp <- prcomp(x[,get_canalysis_variables(s1)], scale=scale, center=center)
   vars <- pcomp[["sdev"]]^2
   vars <- vars/sum(vars)
   first_comps <- (cumsum(vars) < cumproportion)
   pcomp <- predict(pcomp)[,first_comps]
   df2 <- cbind(x, pcomp)
   # SETTINGS
   s2 <- generate_cdata_settings(df2)
   s2[row.names(s1),] <- s1[,colnames(s2)]
   s2[,"in_canalysis"] <- "FALSE"
   new_names <- row.names(na.omit(s2[s2[,"visu_ycoord"] == "NA",]))
   s2[new_names,"in_canalysis"] <- "TRUE"
   s2[new_names,"group"] <- new_names
   s2[new_names,"visu_groups"] <- "PCA"
   s2[new_names,"visu_ycoord"] <- 10*(1:length(new_names)/length(new_names))
   s2[new_names,"fun_transform"] <- "transform_AVG transform_SIGMA"
   max_heatmap <- max(as.numeric(s2[,"heatmap_ycoord"]),na.rm=TRUE)
   s2[new_names,"heatmap_ycoord"] <- as.character((max_heatmap+1):(max_heatmap+length(new_names)))
   x2 <- set_cdata(data=df2, settings=s2, prefix=prefix, data_o=attr(x,'data_o'))
   return(x2)
}

`get_cdata_princomp` <- function(x, cumproportion=0.95, scale=TRUE, center=FALSE){
   # RETRIEVE PREVIOUS PARAMETERS
   s1 <- as.matrix(attr(x,"settings"))
   prefix <- gsub('[0-9]{4}-[0-9]{2}-[0-9]{2}_','',paste2(attr(x,"prefix"),"_ON_PRCOMP_"))
   # PRINCIPAL COMPONENTS, 95%
   pcomp_m <- princomp(x[,get_canalysis_variables(s1)], scale=scale, center=center)
   vars <- pcomp_m[["sdev"]]^2
   vars <- vars/sum(vars)
   first_comps <- names(vars)[cumsum(vars) < cumproportion]
   pcomp <- pcomp_m[["scores"]][,first_comps]
   df2 <- cbind(x, pcomp)
   # SETTINGS
   s2 <- generate_cdata_settings(df2)
   s2[row.names(s1),] <- s1[,colnames(s2)]
   s2[,"in_canalysis"] <- "FALSE"
   new_names <- row.names(na.omit(s2[s2[,"visu_ycoord"] == "NA",]))
   s2[new_names,"in_canalysis"] <- "TRUE"
   s2[new_names,"group"] <- new_names
   s2[new_names,"visu_groups"] <- "PCA"
   s2[new_names,"visu_ycoord"] <- 10*(1:length(new_names)/length(new_names))
   s2[new_names,"fun_transform"] <- "transform_AVG transform_SIGMA"
   max_heatmap <- max(as.numeric(s2[,"heatmap_ycoord"]),na.rm=TRUE)
   s2[new_names,"heatmap_ycoord"] <- as.character((max_heatmap+1):(max_heatmap+length(new_names)))
   x2 <- set_cdata(data=df2, settings=s2, prefix=prefix, data_o=attr(x,'data_o'))
   attr(x2,'princomp') <- pcomp_m
   return(x2)
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
   which_var         <- which(settings[,"in_canalysis"] == "TRUE" | settings[,"in_canalysis"] == " TRUE")
   return(row.names(settings)[which_var])
}

`get_fun_stats` <-
function(fun_name="oddratios",...){
   # CHEMO-INF: CHI2TEST STATS AND RESIDUALS 
   if(fun_name == "chi2test")
      return(function(data,class) { return(stats_chi2test(data,class,...)) })
   # CHEMO-INF: JOINT DISTRIBUTION
   if(fun_name == "jointdistrib")
      return(function(data,class) { return(stats_jointdistrib(data,class,...)) })
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
function(x){
   # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
   bictable     <- as.data.frame(print(attr(x,"bicanalysis"))[["long"]][
                        , c("modelName","G","rseed","relative_BIC")])
   best_bics    <- bictable[na.omit(order(bictable[,"relative_BIC"]
                                ,decreasing=FALSE)[1:attr(x,"nbr_top_models")]),]
   best_models  <- apply(best_bics[,c("modelName","G","rseed")],1,paste,collapse=",",sep="")
   return(best_models)
}

`init_data_cc` <-
function(data,settings)
{
   var_names         <- get_canalysis_variables(settings)
   data_mat          <- data[row.names(na.omit(data[,var_names])),]
   return(data_mat)
}


`get_long2wide_table`          <- function(x, fun_aggregate=mean, 
   var_aggregate="BIC", var_x="G", var_y="modelName"){
   # x <- na.omit(x)
   tmp <- aggregate(x[, var_aggregate], by=list(var_x=x[,var_x],
            var_y=x[,var_y]), FUN=fun_aggregate)
   tmp <- reshape(tmp, idvar="var_x", timevar="var_y", direction="wide")
   row.names(tmp) <- tmp[,"var_x"]
   colnames(tmp) <- gsub(paste2("x."),"",colnames(tmp))
   return(tmp[,-1])
}

`set_bicanalysis` <- function(x){
   bic_out <- bic_rseed <- ranking <- list()
   fun_pattern <- attr(x,"fun_bic_pattern")
   bic_table <- cbind(attr(x, "cfun_params"),BIC=NA,relative_BIC=NA)
   # BIC_TABLE (LONG FORMAT): ABSOLUTE AND RELATIVE
   for(i in 1:length(x))
      bic_table[i,"BIC"] <- attr(x[[i]],"model")[["BIC"]]
   bic_table[,"relative_BIC"] <- 100*(bic_table[,"BIC"]/max(bic_table[,"BIC"],na.rm=TRUE)-1)
   bOut <- structure(as.matrix(bic_table), BIC_stats=NULL, BIC_array=NULL, class="bicanalysis")
   rseed_vector <- unique(bic_table[,"rseed"])
   # COMPUTE BIC STATS PATTERNS ACCORDING TO FUN_PATTERN: MEAN, SD, MEDIAN, ETC.
   if(!is.null(fun_pattern) && length(rseed_vector)>1){
      for(i in names(fun_pattern))
         for(which_bic in list("BIC","relative_BIC"))
            bic_out[[paste2(which_bic, " ", i)]] <- get_long2wide_table(bic_table, 
               var_aggregate=which_bic, fun_aggregate=fun_pattern[[i]])
      # PREPARE AN ARRAY OF BIC TABLES 
      for(i in 1:length(rseed_vector)){
         tmp <- get_long2wide_table(bic_table[bic_table[,"rseed"]==rseed_vector[i],], 
                   var_aggregate="BIC")
         if(i == 1)
            bic_array <- tmp
         else
            bic_array <- abind(bic_array,tmp,along=3)
      }
      dimnames(bic_array)[[3]]  <- rseed_vector
      # GIVEN A modelName AND G, AVERAGE RANK OF EACH RSEED
      ranking[["Average rank of each rseed"]] <- sort(apply(apply(abs(bic_array),1:2,rank),1,patternMean))
      # GIVEN A modelName, RANK OF G ACCROSS ALL THE RANDOM SEEDS
      ranking[["Given a modelName, average rank of G"]] <- apply(apply(abs(bic_array),2:3,rank),1:2,patternMean)
      # GIVEN A G, RANK OF modelName ACCROSS ALL THE RANDOM SEEDS
      tmp_array <- apply(abs(bic_array[,,1]),1,rank)
      for(idx in 2:dim(bic_array)[[3]])
         tmp_array <- abind(tmp_array, apply(abs(bic_array[,,idx]),1,rank),along=3)
      ranking[["Given a G, average rank of modelName"]] <- apply(tmp_array,1:2,patternMean)
      ranking[["In all rseeds, the best one"]] <- matrix(NA, nrow(bic_array), ncol(bic_array), 
                                   dimnames=list(dimnames(bic_array)[[1]], dimnames(bic_array)[[2]]))
      ranking[["In all rseeds, the best BIC"]] <- ranking[["In all rseeds, the best one"]]
      for(i in 1:nrow(ranking[["In all rseeds, the best one"]])){
         for(j in 1:ncol(ranking[["In all rseeds, the best one"]])){
            best_seed <- best_score <- NA
            if(length(na.omit(bic_array[i,j,]))>0){
               best_seed <- dimnames(bic_array)[[3]][which(bic_array[i,j,] == patternMax(bic_array[i,j,]))]
               best_score<- patternMax(bic_array[i,j,])
            }
            ranking[["In all rseeds, the best one"]][i,j] <- best_seed
            ranking[["In all rseeds, the best BIC"]][i,j] <- best_score
         }
      }
      attr(bOut,'BIC_stats') <- append(bic_out,ranking) 
      attr(bOut,'BIC_array') <- bic_array
   }
   if((!is.null(fun_pattern)) && length(rseed_vector)==1){
      ranking[["Given a G, average rank of modelName"]] <- apply(abs(bic_out[["BIC mean"]]),1,rank)
      ranking[["Given a modelName, average rank of G"]] <- apply(abs(bic_out[["BIC mean"]]),2,rank)
      attr(bOut,'BIC_stats') <- ranking 
   }
   attr(x,"bicanalysis") <- bOut
   return(x)
}

`print.bicanalysis` <- function(x,silent=TRUE,...){
   x_out <- list()
   x_out[["long"]] <- x[1:nrow(x),]
   for(a_n in names(attr(x,"BIC_stats")))
      x_out[[a_n]] <- attr(x,"BIC_stats")[[a_n]]
   x_out[["BIC_array"]] <- attr(x,"BIC_array")
   if(!silent)
      print(as.data.frame(x[1:nrow(x),]))
   return(x_out)
}

`initPlot` <- 
function(cdata, figIdx, plotNbr=1, type='MixtModel'){
   # CREATE FIGURE DIRECTORY
   dir.create(attr(cdata,'figDir'),showWarnings=FALSE)
   figFile <- paste2(attr(cdata,'figDir'),attr(cdata,"prefix"),'_',sprintf('%03d',figIdx),'-',type,'.pdf')
   print(figFile)
   # NUMBER OF SUBPLOT IN A FIGURE
   plotI <- plotJ <- as.integer(sqrt(plotNbr))
   if(plotI^2 < plotNbr)
      plotJ <- plotJ+1 
   # FIGURE OUTPUT
   maiParam <- attr(cdata,'mai')
   if(type=='BB')
      maiParam <- maiParam+c(0.5,0,0.5,0) 
   pdf(figFile)
   par(mfrow = c(plotI,plotJ), mai=maiParam)
   figIdx <- figIdx+1
   return(figIdx)
}

`closePlot` <-
function(){
   dev.off()
}

`plot.cdata` <- 
function(x, ...) {
   xSettings <- as.matrix(na.omit(attr(x,"settings")[,c("visu_groups","visu_ycoord")]))
   splitGroup <- split(row.names(xSettings),xSettings[,"visu_groups"])
   figIdxBB <- figIdxH <- 1
   for(i in 1:length(splitGroup)){
      # PRODUCE BOXPLOT BY SPLIT GROUP
      figIdxBB <- initPlot(x, figIdxBB, plotNbr=1, type='BB')
      boxplot(as.data.frame(x[,splitGroup[[i]]]),'las'=2, 
         main=paste2(" Boxplot ",names(splitGroup)[i]) , cex.axis=1)
      closePlot()
      figIdxH <- initPlot(x, figIdxH, plotNbr=length(splitGroup[[i]]), type='H')
      for(varName in splitGroup[[i]]){
         # PREPARE THE TEXT FOR EACH HISTOGRAM
         vIdx <- which(names(attr(x,"tdata")) == varName)
         xlabText <- paste2(varName,' ', length(table(x[,varName])), " values\n")
         if(length(vIdx) == 1){
            for(tn in names(attr(x,"tdata")[[vIdx]])){
               if(tn == "transform_adjust")
                  xlabText <-  paste2(xlabText, tolower(gsub("transform_","",tn))
                     , "=", attr(x,"tdata")[[vIdx]][[tn]][["print"]], " ")
               else
                  xlabText <-  paste2(xlabText, tolower(gsub("transform_","",tn))
                     , "=", sprintf("%.2f",attr(x,"tdata")[[vIdx]][[tn]])," ")
               }
            }
         # PRODUCE HISTOGRAM BY SPLIT GROUP
         hist(x[,varName], main=NULL,xlab=xlabText, ylab="", cex=0.5, las=2,'new'=TRUE, border=NULL)
      }
      closePlot()
   }
}

`plot.cresult` <- 
function(x, query = NULL, ...) {
   if(is.null(query))   
      q <- 1:length(x)
   else 
      q <- query
   nbrPlot <- length(attr(x,"fun_plot"))
   for(i in q){
      if(!is.na(attr(x[[i]],"model")$loglik)){
         par(mfrow=mfrow, 'las'=1,'mai'=c(0.2,0.7,0.5,0.7),new=FALSE)
         for(fun_plot in unlist(attr(x,"fun_plot")))
            try(fun_plot(x[[i]]))
      }
   }
}

`html_format_cells` <- function(x, fmt="%.2f", quant_threshold=quantile(c(-1,1),probs=seq(0,1,0.333))){
   x_nbr <- as.numeric(x) 
   col_vector <- brewer.pal(length(quant_threshold)+1,"OrRd")
   if(!is.na(x_nbr)){
      #x <- gsub("\\.00","",gsub("^0\\.",".",sprintf(fmt,x)))
      x <- sprintf(fmt,x)
      for(i in 1:length(quant_threshold)){
         if(x_nbr <= quant_threshold[i])
            x <- paste2("<font style='background-color:",col_vector[1+i],";text-align:right;'>",x,"</font>")
         if((i == length(quant_threshold)) && (x_nbr >= quant_threshold[i]))
            x <- paste2("<font style='background-color:",col_vector[1+i],";text-align:right;'>",x,"</font>")
      }
   }
   else
      x <- ""
   return(x)
}

`html_format_numeric` <- function(x, fmt="%.2f", col_threshold=c(-1,1)){
   x_nbr <- as.numeric(x) 
   if(!is.na(x_nbr)){
      x <- paste2("<div align='right'>" 
                , gsub("\\.00","",gsub("^0\\.",".",sprintf(fmt,x))) 
                , "</div>")
      if(x_nbr >= col_threshold[2])
         x <- paste2("<font color='red'>",x,"</font>")
      if(x_nbr <= col_threshold[1])
         x <- paste2("<font color='blue'>",x,"</font>")
   }
   else
      x <- ""
   return(x)
}

`print.cresult` <- 
function(x, html = TRUE, img=FALSE, img_size=c(1024,1024),...) {
   x_comp <- compare_cresult(x)
   x_out <- print(attr(x,"bicanalysis"))
   x_out[["Best models"]] <- get_best_models(x)
   x_csv <- x_out 
   if(html){
      for(p_name in names(x_comp)){
         pxcomp <- as.matrix(x_comp[[p_name]][["full"]])
         pxcomp2 <- as.matrix(x_comp[[p_name]][["html"]])
         row.names(pxcomp) <-row.names(pxcomp2) <- gsub("_[-]?0\\.0","",row.names(pxcomp2))
         colnames(pxcomp) <- colnames(pxcomp2) <- gsub("_[-]?0\\.0","",colnames(pxcomp2))
         x_out[[p_name]] <- pxcomp2
         x_csv[[p_name]] <- pxcomp
      }
      f_out <- paste2(attr(x,"prefix"),"_report")
      HTMLStart(outdir=getwd(), filename=f_out, HTMLframe=TRUE
         , CSSFile="R2HTML.css", Title=attr(x,"prefix"))
      for(i in sort(names(x_out))){
         a_name <- gsub(" ","%20",i)
         HTML.title(paste2("<a name='",a_name,"'>&nbsp;</a>",i))
         file_csv <- paste2(attr(x,"prefix"),"_",i,".csv")
         HTML(paste2("<a href=",f_out,"_main.html#",a_name," target=main>",i,"</a>", 
            " (<a href='./",file_csv,"' target=main>CSV file</a>)"), file=paste2(f_out,"_menu.html"))
         write.table(x_csv[[i]], na="", dec=",", sep=";", file=file_csv)
         HTML(x_out[[i]],align="left")
      }
      if(img){
         for(m_name in names(x)){
            file_img <- plot(x,query = m_name,device="PNG", title=m_name)
            HTMLInsertGraph(file_img,Caption = paste2("Figure. Graphic characteristics of model ",m_name)
                  , GraphBorder = 0, WidthHTML=img_size[1],HeightHTML=img_size[2])
         }
      }
      HTMLStop()
   }
   else
      return(x_csv)
}

`fun_cmodel` <-
function(x) {
   # RETRIEVE DATA FOR THE CANALYSIS
   cdata        <- get_cdata(attr(x,"cdata"))
   for(i in 1:length(x)){
      cat("\n",names(x)[[i]])
      attr(x[[i]],"model") <- x[[i]](x[[i]],cdata)
      if(!is.na(attr(x[[i]],"model")$loglik)){
         cdata_l <- attr(x,"cdata")[,order(attr(attr(x,"cdata"),"settings")[,"heatmap_ycoord"],na.last=NA)] 
         cdata_l <- cbind(cdata_l,class=attr(x[[i]],"model")[["labelling"]])
         cat("-> Patterns")
         # COMPUTE PATTERNS
         attr(x[[i]],"pattern") <- NULL
         for(p in names(attr(x,"fun_pattern")))
            attr(x[[i]],"pattern")[[p]] <- compute_pattern(cdata_l,attr(x,"fun_pattern")[[p]])
         # FOR ORDERING, MAKE A DENDROGRAM FROM THE CENTER PATTERN (NUMBER 1)
         if(!is.null(attr(x[[i]],"pattern"))){
            cat("-> Dendros")
            avg_p <- attr(x[[i]],"pattern")[[1]]
            attr(x[[i]],"dendro_cluster") <- hclust(dist(avg_p))
            attr(x[[i]],"dendro_var") <- hclust(dist(t(avg_p)))
            # REORDER THE PATTERNS
            cat("-> Ordering")
            for(p in names(attr(x,"fun_pattern")))
               attr(x[[i]],"pattern")[[p]] <- attr(x[[i]],"pattern")[[p]][attr(x[[i]],"dendro_cluster")$order,]
   #                                                , rev(sort(colnames(avg_p)))]
            attr(x[[i]],"pattern")[["class_count"]] <- table(cdata_l[,"class"])[attr(x[[i]],"dendro_cluster")$order]
            # SELECT DIVERGING COLORS FROM RColorBrewer
            attr(x[[i]],"cluster_colors") <- brewer.pal(nrow(avg_p),"Set1")
         }
         # STATS
         fun_stats <- attr(x,"fun_stats")
         cat("-> Stats")
         for(n in names(fun_stats))
            attr(x[[i]],"stats")[[n]] <- fun_stats[[n]](data=attr(x,"cdata"), class=attr(x[[i]],"model")[["z"]])
      }
   }
   # WRITE A BIC SCORE TABLE AS AN OUTPUT
   x <- set_bicanalysis(x)
   return(x)
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
      csv_output        <- paste2(attr(x,"prefix"), '_'
                                , gsub(",","_",q),".csv") 
      print(csv_output)
      write.csv2(out,file=csv_output)
   }
}

`patternMean` <- function(x){
   return(mean(x,na.rm=TRUE))
}

`patternMedian` <- function(x){
   return(median(x,na.rm=TRUE))
}

`patternLowquant` <- function(x){
   return(quantile(x,probs = 0.025,na.rm=TRUE))
}

`patternUpquant` <- function(x){
   return(quantile(x,probs = 0.975,na.rm=TRUE))
}

`patternMax` <- function(x){
   return(max(x,na.rm=TRUE))
}

`patternSd` <- function(x){
   return(sd(x,na.rm=TRUE))
}

`set_cresult` <-
function(cdata=NULL, cfun=fun_mbc_em, cfun_settings=list(modelName=c("EII",
   "VII"), G=3:5, rseed=6013:6015), fun_pattern=list(mean=patternMean,
   median=patternMedian, lowquant=patternLowquant, upquant=patternUpquant) ,
   fun_plot=list(plot_parcoord=get_plot_fun(type="plot_parcoord"),
   plot_legend=get_plot_fun(type="plot_legend"),
   plot_image=get_plot_fun(type="plot_image"),
   plot_dendro_cluster=get_plot_fun(type = "plot_dendro_cluster"),
   plot_dendro_var=get_plot_fun(type="plot_dendro_var")),
   fun_stats=list(oddratios=get_fun_stats(fun_name="oddratios")),
   fun_bic_pattern=list(mean=patternMean, median=patternMedian,
   lowquant=patternLowquant, upquant=patternUpquant, sd=patternSd),
   nbr_top_models=5){
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
      fun_call          <- function(cmodel,cdata){return(do.call(cfun,append(list(data=cdata),attr(cmodel,"cfun_settings"))))}
      cresult_out[[id_name]] <- structure(fun_call,
         cfun_settings=param_set_df[i,], model=NULL, pattern=list(),
         dendro_cluster=NULL, dendro_var=NULL, cluster_colors=NULL,
         stats=list(), class="cmodel")
   }
   # PRE-EVAL FUN_PLOT TO RETRIEVE THE DIFFERENT PLOTTING FUNCTIONS
   fun_plot2 <- list()
   for(lfun in unlist(fun_plot))
      fun_plot2 <- c(fun_plot2, lfun(cdata))
   # CONSTRUCT A 'CRESULT' DATA STRUCTURE WITH CMODELS IN IT
   cresult_out <- structure(cresult_out, cdata=cdata,
      prefix=attr(cdata,"prefix"), cfun_params=param_set_df,
      fun_plot=fun_plot2, fun_stats=fun_stats, fun_pattern=fun_pattern,
      fun_bic_pattern=fun_bic_pattern, nbr_top_models=nbr_top_models,
      bicanalysis=NULL, rinfo=sessionInfo(), class = "cresult")
   return(cresult_out)
}

`set_cdata` <-
function(data, data_o = NULL, tdata = NULL, settings, init_fun = list(init_data_cc), prefix =""){
   tdata_list <- list()
   settings <- settings[row.names(settings)[!is.na(settings[,"group"])],]
   # BACKUP DATA INTO DATA_O 
   if(is.null(data_o) )
      data_o <- data
   # DO INIT FUNCTIONS OF THE DATA
   for(fun in init_fun)
      data <- fun(data,settings)
   # TRANSFORM THE DATA
   data <- data[,row.names(settings)]
   t_result <- transform_cdata(data,settings,tdata)
   data_out <- data.matrix(t_result[["data"]])
   tdata_out <- t_result[["tdata"]]
   cdata_out <- structure(data_out
                        , data_o        = data_o  
                        , tdata         = tdata_out
                        , settings      = settings 
                        , prefix        = paste2(today(),"_",prefix)
                        , xlim          = c(-3,3)
                        , ylim          = c(0,50)
                        , mai           = c(0.6,0.3,0.05,0.05)
                        , figDir  = 'figures/' 
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
      T_s <- as.matrix(T_s)
      min2 <- function(x){if(min(x)<0){return(min(x))}else{return(0)}}
      max2 <- function(x){if(max(x)>10){return(max(x))}else{return(10)}}
      ylim<- c(min2(as.numeric(T_s[,"visu_ycoord"])),max2(as.numeric(T_s[,"visu_ycoord"])))
      plot(x=0 ,new=TRUE ,ann=FALSE ,pch=18,col="white",axes=FALSE ,xlim=xlim , ylim=ylim)
      title(main = paste2("(",paste(as.matrix(attr(cmodel, "cfun_settings")),collapse=","),") ",title), cex = cex * 0.8)
      for(x_name in pattern){
         pattern <- attr(cmodel,"pattern")[[x_name]][,row.names(T_s)]
         lwd <- 3
         lty <- "solid"
         if(x_name != "mean") {
            lwd <- 1
            lty <- "dashed"
         }
         for(g in 1:nrow(pattern)){
            D_i_s <- cbind(X=as.numeric(pattern[g,]),Y=as.numeric(T_s[,"visu_ycoord"]))
            D_i_s <- D_i_s[sort.list(D_i_s[,"Y"]),]
            gap <- 1.05*abs(median(D_i_s[1:(nrow(D_i_s)-1),"Y"] - D_i_s[2:nrow(D_i_s),"Y"]))
            for(l in 1:(nrow(D_i_s)-1))
               if(!(D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap))
                  arrows(D_i_s[l,1], D_i_s[l,2], D_i_s[l+1,1], D_i_s[l+1,2], col=attr(cmodel,"cluster_colors")[g], length=0, lwd=lwd, lty=lty)
               # ELSE, AS (is_white_gap == TRUE) THEN DO NOT DRAW ANY ARROW...
         }
      }
      fv_idx <- 1:nrow(T_s)
      if(nrow(T_s) > 20)
         fv_idx <- as.integer(quantile(1:nrow(T_s)))
      axis(2,at=as.numeric(T_s[fv_idx,"visu_ycoord"]),labels=colnames(pattern)[fv_idx],las=2,tick=FALSE)
      axis(1,at=seq(from = xlim[1] , to = xlim[2], by = (xlim[2] - xlim[1])/4))
   })
} 

`get_plot_fun` <- function(type="plot_parcoord", title=NULL, xlim=c(-3, 3), zlim=c(-2,2),
   xy=c(-2.2, 0), cex=0.9, pattern="mean",
   color_gradient=rev(brewer.pal(9,"RdBu")),range_fv=NULL){
#   color_gradient=colorRampPalette(rev(brewer.pal(11,"RdBu")))(8)){
   return(function(cdata){
      return(get_plot_fun2(type=type, cdata=cdata, xlim=xlim, xy=xy, cex=cex,
         title=title, pattern=pattern, color_gradient=color_gradient,range_fv=range_fv, zlim=zlim))
      })
}

`get_plot_fun2` <- function(type, cdata, title, xlim, xy, cex, pattern,
   color_gradient,range_fv,zlim){
   if(type == "plot_parcoord"){
      cfun_parcoord     <- c()
      for(v in unique(as.character(attr(cdata,"settings")[,"visu_groups"]))){
         if(!is.na(v)){
            T_s            <- na.omit(attr(cdata,"settings")[attr(cdata,"settings")[,"visu_groups"] == v,])
            lfun           <- eval(call("plot_parcoord",xlim = xlim , cex = cex , title = v, T_s = T_s , pattern = pattern))
            cfun_parcoord  <- c(cfun_parcoord, lfun)
         }
      }
      return(cfun_parcoord)
   }
   if(type == "plot_series"){ 
      if(is.null(range_fv))
         range_fv <- 1:ncol(cdata)
      range_fv <- colnames(cdata)[range_fv]
      return(
         function(cmodel, apattern = pattern, acex = cex, axlim=xlim,rangeFeatVector=range_fv){
               x <- attr(cmodel,"pattern")[[apattern]][,rangeFeatVector]
               plot(ts(t(x)), plot.type='single', col=attr(cmodel,"cluster_colors"), lwd=2, axes=FALSE, 
                  ylim=axlim, xlab='',ylab='', main=paste(as.matrix(attr(cmodel, "cfun_settings")),sep=',',collapse=','))
               fv_idx <- 1:ncol(x)
               if(ncol(x) > 20)
                  fv_idx <- as.integer(quantile(fv_idx))
               axis(1, at=(1L:ncol(x))[fv_idx], labels = colnames(x)[fv_idx], las = 2, line = -0.5, tick = TRUE, cex.axis = acex)
               axis(2, at=as.integer(seq(axlim[1],axlim[2],length.out=4)),las=2)
               d <- table(attr(cmodel,'model')$labelling)[row.names(x)]
               d <- apply(cbind(names(d),' (',d,')'),1,paste2)
               legend("bottomleft",legend=d,col=attr(cmodel,"cluster_colors"),lwd=2, bty='n')
            })
            }
   if(type == "plot_image")            
      return(
         function(cmodel, apattern = pattern, acolgrad = color_gradient, atitle = title, acex = cex, azlim=zlim){
               x <- attr(cmodel,"pattern")[[apattern]]
               image(1L:nrow(x), 1L:ncol(x), x, xlim = 0.5 + c(0, nrow(x)),
                  ylim = 0.5 + c(0, ncol(x)), axes = FALSE, xlab = "", ylab = "",
                  col=acolgrad, main=atitle,add=FALSE, zlim=azlim)
               fv_idx <- 1:ncol(x)
               if(ncol(x) > 20)
                  fv_idx <- as.integer(quantile(fv_idx))
               axis(2, at=(1L:ncol(x))[fv_idx], labels = colnames(x)[fv_idx], las = 2, line = -0.5, tick = 0, cex.axis = acex)
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
                     if(type == "plot_dendro_var") 
                        x <- attr(cmodel,"dendro_var")
                     else
                        x <- attr(cmodel,"dendro_cluster")
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
   b <- as.data.frame(cbind(data,class = class))
   b <- reshape(b[,c("family","member","class")],timevar="member",idvar="family",direction="wide")
   b <- b[,-1]
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
   s <- t(s[,1:length(class_set)])
   s_out <- apply(s, 1:2, html_format_numeric, fmt='%.1f')
   return(list(out=s_out,s))
}

`stats_logodds` <-
function(data,class,fun_midthreshold=median){
   #
   nclass <- ncol(class)
   class <- map(class, warn=FALSE)
   ldata <- cbind(data,class=class)
   sgroup <- as.data.frame(na.omit(attr(data,"settings")[,c("visu_groups","group")]))
   s <- matrix(0,nclass, length(unique(sgroup[,"group"])),
                dimnames=list(1:nclass, sort(unique(sgroup[,"group"]))))
   class_vector <- 1:nclass
   gr_vector <- sort(unique(as.character(sgroup[,"group"])))
   nrow2 <- function(x){y <- nrow(x) ; if(is.null(y)){return(0)}else{return(y)}}
   sum2 <- function(x){if(!is.null(dim(x))){return(apply(x,1,sum))}else{return(x)}}
   # FOR EACH HIERARCHICAL SUBSET OF OUTCOMES, MAKE SUM SCORES
   for(gr in gr_vector){
      m <- cbind(SScore=sum2(ldata[,row.names(sgroup[sgroup[,"group"]==gr,])]),class=class)
      # FOR EACH PATTERN, DERIVE ITS STATISTIC
      for(g in class_vector){
         middleScore <- fun_midthreshold(m[,"SScore"],na.rm=TRUE)
         m11 <- nrow2(m[m[,"class"]==g & m[,"SScore"] >= middleScore,])
         m12 <- nrow2(m[m[,"class"]==g & m[,"SScore"] < middleScore,])
         m21 <- nrow2(m[m[,"class"]!=g & m[,"SScore"] >= middleScore,])
         m22 <- nrow2(m[m[,"class"]!=g & m[,"SScore"] < middleScore,])
         mat <- matrix(c(m11,m12,m21,m22),2,2,byrow=TRUE)
         # LOG OF THE ODD-RATIO ALSO NAMED CROSS PRODUCT
         s[g,match(gr,gr_vector)] <- log(mat[1,1] * mat[2,2] / (mat[1,2] * mat[2,1]))
         # dimnames(s)[[2]][match(gr,gr_vector)] <- paste(gr,sprintf("%.1f",middleScore),sep="_",collapse="")
      }
   }
   s_out <- apply(s, 1:2, html_format_numeric, fmt='%.1f')
   return(list(out=s_out,s))
}

`stats_chi2test` <- 
function(data, class, class_col='Class', probs=seq(0,1,0.4)[2:3]){
   ldata <- cbind(data, modelClass=map(class))
   ldata <- cbind(ldata, target=as.character(attr(data,'data_o')[row.names(ldata),class_col]))
   s_test <- chisq.test(xtabs( ~.,ldata[,c('modelClass',"target")]), simulate.p.value=TRUE)
   s <- cbind(s_test[["residuals"]]^2, chi2test=NA)
   s[1:2,"chi2test"] <- c(sprintf('%.1e',s_test[["p.value"]]), sprintf('%.1f',s_test[["statistic"]]))
   colnames(s)[match("chi2test",colnames(s))] <- paste2(class_col, " (residual)")
   s_out <- apply(s, 1:2, html_format_cells, fmt='%.1f' 
      , quant_threshold=quantile(as.numeric(s_test[["residuals"]]^2), probs=probs))
   return(list(out=s_out,s))
}

`stats_jointdistrib` <- 
function(data, class, class_col='Class'){
   ldata <- cbind(data, modelClass=map(class))
   ldata <- cbind(ldata, target=as.character(attr(data,'data_o')[row.names(ldata),class_col]))
   s <- xtabs( ~.,ldata[,c("modelClass","target")])
   s_out <- apply(s, 1:2, html_format_numeric, fmt='%d')
   return(list(out=s_out,s))
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
   tfun         <- list(function(x){sd_est <- sd(x,na.rm=TRUE) ; if(sd_est>0) return(sd_est) else return(1)}
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
   # BEFORE TRANSFORMING THE DATA, CONVERT IT TO NUMERIC
   data <- data.matrix(data)
   if(is.null(tdata)) tdata <- list()
   settings <- as.matrix(settings)
   # FOR EACH VAR
   for(var in row.names(settings)){
      tdata_idx                 <- which(names(tdata) == var)
      if(length(tdata_idx) == 0) 
         tdata[[var]]  <- list()
      # USE (MULTIPLE)-SPACES AS SEPARATOR BETWEEN TRANSFORMATION FUNCTIONS
      fun_transform             <- strsplit(settings[var,"fun_transform"],"[ ]+")[[1]]
      for(lfun in fun_transform){
         if(is.na(lfun)) 
            break ;
         # IN CASE 'LFUN' AS AN INSIDE PARAMETER, RETRIEVE IT (1) REPLACE ALL
         # ('\") BY SPACES AND THEN (2) SPLIT
         lfun                   <- gsub(" ","",lfun,extended=FALSE)
         lfun                   <- gsub("\"","'",lfun,extended=FALSE)
         lfun                   <- strsplit(gsub("')","",gsub("('"," ",lfun,extended=FALSE),extended=FALSE),"[ ]+")[[1]]
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

`perform_chi2_tests` <- 
function(m,tonominal=FALSE,simulate.p.value = TRUE, B = 10000){
   m2 <- m
   if(tonominal)
      m2 <- cbind(m,sup_med=as.numeric(m[,1])>median(as.numeric(m[,1])))[,3:2]
   a <- xtabs(formula=~.,data=m2)
   if(nrow(a) > 2) range_row <- 1:nrow(a)
   else range_row <- 1
   apvalues <- matrix(NA,max(range_row),ncol(a))
   for(i in range_row){
      for(j in 1:ncol(a)){
         a2by2 <- matrix(c(a[i,j],sum(a[i,][-j]),sum(a[,j][-i]),sum(a[-i,-j])),2,2,byrow=TRUE)
         # FOR TESTING: cat("\n",i,j) print(a2by2)
         chi2_test <- chisq.test(a2by2,simulate.p.value = simulate.p.value, B = B)
         apvalues[i,j] <- sign(chi2_test$residuals[1,1]) * chi2_test$p.value
      }
   }
   if(nrow(apvalues) > 1)
      row.names(apvalues) <- row.names(a)
   else if(nrow(apvalues) == 1 && tonominal) 
      row.names(apvalues)[1] <- paste2('> med ',colnames(m)[1])
   else 
      row.names(apvalues)[1] <- row.names(a)[1]
   return(apvalues)
}

`word_cloud_analysis` <- 
function(m, convert2nominal=NA,title=NA,filename="WC",test_threshold=0.05,magn=0.05,css_style="font-family:times;",colors=c('#fddbc7','#d1e5f0')){
   HTMLStart(outdir="./",filename=paste2(today(),filename),HTMLframe=FALSE, Title=title,autobrowse=FALSE)
   HTML(title)
   a <- htmlTable <- matrix("",0,length(unique(m[,"class"])))
   for(f in colnames(m[,-grep("class", colnames(m))])){
      tonominal <- f %in% convert2nominal
      X2 <- perform_chi2_tests(m[,c(f,"class")], tonominal=tonominal)
      a <- rbind(a, X2)
      line <- matrix("",1,ncol(htmlTable))
      for(j in 1:ncol(X2)){
         for(i in 1:nrow(X2)){
            if(abs(X2[i,j]) <= test_threshold){
               line[j] <- paste2(line[j],"<font style='",css_style,
                                 ";background-color:",colors[as.integer(1.5+sign(X2[i,j])/2)],
                                 ";font-size:", -log(abs(X2[i,j]))*magn, "cm'>", 
                                 row.names(X2)[i], "</font><br>")
            }
         }
      }
      htmlTable <- rbind(htmlTable,line)
   }
   a <- matrix(as.numeric(a),nrow(a),ncol(a),dimnames=list(sub("> med ","",row.names(a)),colnames(a)))
   row.names(htmlTable) <- colnames(m[,-grep("class", colnames(m))])
   colnames(htmlTable) <- 1:ncol(htmlTable)
   HTML(htmlTable)
   HTMLStop()
   return(a)
}

`plot.association_analysis` <- 
function(x, range_fv=NULL, p_quantile=0.95, ylim=NA, ylim2=c(15,0), xlab="", text="", ...){
   postscript(paste2(attr(x,'prefix'),'.ps'))
   par('mfrow'=c(3,2),mai=c('0.5','0.5','0.5','1'))
   for(i in 1:length(attr(x,'fv'))){
      if(length(range_fv)>0)
         a <- cbind(scores=attr(x,'sc_med')[[i]][range_fv], pvalues=apply(attr(x,'x2')[[i]][,range_fv],2,quantile,probs=p_quantile))
      else
         a <- cbind(scores=attr(x,'sc_med')[[i]], pvalues=apply(attr(x,'x2')[[i]],2,quantile,probs=p_quantile))
      if(length(ylim) != 2)
         ylim <- as.integer(range(a[,'scores']))
      title <- paste2(text, ' (#',i,')')
      plot_patternAndPvalues(a, title=title, ylim2=ylim2, ylim=ylim, xlab=xlab)
   }
   graphics.off()
}


`plot_patternAndPvalues` <-
function(a, title="", ylim=c(-2,38), ylim2=c(15,0), xlab=""){
   # P-VALUES
   plot(x=0, ann=FALSE, pch=18, col="white", axes=FALSE, xlim=c(1,nrow(a)), ylim=ylim2)
   lines(list(x=1:nrow(a),y=-log(a[,'pvalues'])), col='grey', lty=1, lwd=1)
   lines(list(x=c(5,nrow(a)), y=-log(c(0.05,0.05))),col='grey',lwd=1)
   lines(list(x=c(5,nrow(a)), y=-log(c(0.01,0.01))),col='grey',lwd=1)
   axis(4, at=-log(c(1,0.05,0.01)), las='1', labels=c('1   ',' .05',' .01'))
   # IMAGE
   b <- matrix(0,nrow(a),1)
   b[which(abs(a[,'pvalues'])<0.05),] <- 0.5
   b[which(abs(a[,'pvalues'])<0.01),] <- 1
   image(x=1:nrow(a),y=-log(c(1,0.05)),z=b,col=brewer.pal(5,"Blues"), breaks=c(-0.1,0.1,0.3,0.6,0.8,1),add=TRUE)
   par(new=TRUE)
   plot(x=0, new=TRUE, ann=FALSE, pch=18, col="white", axes=FALSE, xlim=c(1,nrow(a)), ylim=ylim,xlab=xlab)
   # SCORES
   lines(list(x=1:nrow(a), y=a[,'scores']), lwd=2)
   # RANGES OF FREQUENCIES: HIGHLY DISCRIMINATIVES
   topFeats <- as.numeric(which(abs(a[,'pvalues'])<0.05))
   k <- j <- topFeats[1]
   for(i in topFeats){
      if(i>j+1) k <- c(k,j,i)
      j<- i
   }
   if(!is.integer(length(k)/2))
   k <- c(k,k[length(k)])
   rangeFeats <- matrix(k,length(k)/2,2,byrow=TRUE)
   # SCORE AXES  
   axis(1, at=c(1,unique(rangeFeats),nrow(a)), labels=row.names(a)[c(1,unique(rangeFeats),nrow(a))],las='2')
   seqY <- seq(ylim[1],ylim[2],length.out=4)
   axis(2, at=seqY, labels=sprintf('%1.1e',seqY), las=2)
   title(title)
}


`bootstrap_pvalues` <- 
function(x, B_max=1000,rseed=6013){
   set.seed(rseed)
   master_bootstrap <- sample(1:nrow(x),B_max*nrow(x),replace=TRUE)
   x2 <- matrix(NA,B_max,ncol(x)-1,dimnames=list(list(),colnames(x[,-ncol(x)])))
   for(b in 1:B_max){
      idx <- master_bootstrap[(1+(b-1)*(nrow(x)-1)):(b*(nrow(x)-1))]
      cat(".",b)
      for(f in 1:(ncol(x)-1)){
         x2[b,f] <- chisq.test(xtabs(formula=~.,data=x[idx,c(f,ncol(x))]),simulate.p.value=TRUE, B = 1000)$p.value
         }

   }
   return(x2)
}


`to_binary` <- 
function(x){
   x_bin <- x
   for(i in 1:ncol(x))
      x_bin[,i] <- as.numeric(x[,i])>median(as.numeric(x[,i]))
   return(x_bin)
}


`association_analysis` <- 
function(x, q=NA, B=4){
   if(class(x) == "cresult"){
      # RETRIEVE APPROPRIATE CLASS/SUBTYPE LABELS 
      if(!is.na(q) && length(which(q %in% names(x))) == 0)
         l <- as.character(attr(attr(x,'cdata'),'data_o')[row.names(attr(x,'cdata')),q])
      else{
         if(is.na(q))
            q <- get_best_models(x)[1]
         l <- as.character(attr(x[[q]],'model')$labelling)
      }
      m <- data.frame(attr(attr(x,'cdata'),'data_o')[,colnames(attr(x,'cdata'))],class=l)
      p <- attr(x, 'prefix')
   }
   else{
      p <- today()
      m <- x 
      l <- as.character(x[,'class'])
   }
   # INIT DATA
   prefix <- paste2(p,'_discriminative_features_',gsub(",","_",q))
   mBin <- data.frame(to_binary(m[,-ncol(m)]),class=l)
   x2 <- fv <- sc_med <- list()
   # 
   for(i in sort(unique(l))){
      mTmp <- data.frame(mBin[,-ncol(mBin)], class=(l==i))
      x2[[i]] <- bootstrap_pvalues(mTmp, B_max=B)
      fv[[i]] <- names(which(apply(x2[[i]], 2, quantile, probs=0.5)<0.05))
      sc_med[[i]] <- apply(m[(l==i), -which(colnames(m) %in% c('class',q))],2,median)
   }
   obj_out <- structure(fv, x2=x2, fv=fv, sc_med=sc_med, prefix=prefix, class="association_analysis")
   if((class(x) == "cresult") && (length(which(q %in% names(x)))>0)){
      attr(x[[q]],'association_analysis') <- obj_out
      return(x)
   }
   else
      return(obj_out)
}

`bestNFeatures` <- 
function(x, N=10){
   # FEATURE SELECTION BY RANKING OF THE P-VALUES
   cList <- names(attr(x,'x2'))
   m <- attr(x,'x2')[[1]]
   pTable <- rTable <- fSelFreq <- matrix(NA, ncol(m),length(cList), dimnames=list(colnames(m),cList)) 
   maxP <- matrix(NA,length(cList),1,dimnames=list(cList,'maxP'))
   # NFREQ CONTROLS THE QUANTITY OF FEATURE SELECTED BY CLASS
   nPerClass <- as.integer(N/length(cList))
   fList <- NULL
   for(i in cList){
      pTable[,i] <- abs(attr(x,'sc_med')[[i]])
      rTable[,i] <- order(attr(x,'sc_med')[[i]])
      fSelFreq[,i] <- order(attr(x,'sc_med')[[i]])<=nPerClass
      fList <- unique(c(fList,names(which(fSelFreq[,i]))))
   }
   j <- 1
   while(length(fList)<N){
      i <- 1+j%%length(cList)
      fSelFreq[,i] <- order(attr(x,'sc_med')[[i]]) <= (nPerClass+as.integer(1+j/length(cList)))
      fList <- unique(c(fList,names(which(fSelFreq[,i]))))
      j <- j+1
   }
   for(i in cList)
      maxP[i]<- max(pTable[which(fSelFreq[,i]),i])
   attr(x,'bestNFeatures') <- structure(fList, N=N, pTable=pTable, rTable=rTable, maxP=maxP,class='bestNFeatures')
   return(x)
}


`analysisVarFeat` <- 
function(m, a=NULL, maxFeat=NA, postscript=TRUE,clusterModel='VVI', clusterG=6,
   clusterRseed=6013:6063){
   nFeatMax <- min(ncol(m)-1,maxFeat,na.rm=TRUE)
   iterMax <- as.integer((2*(log(nFeatMax))/log(2)-4))
   iterSeq <- as.integer(2^(2+(1:iterMax)/2))
   x <- NULL
   if(is.null(a))
      a <- association_analysis(m)
   for(i in iterSeq){
      cat('\n',i,' features...')
      aTmp <- bestNFeatures(a,N=i)
      # PREPARE SETTINGS
      settings <- generate_cdata_settings(m)
      settings[,'in_canalysis'] <- FALSE 
      settings[as.character(attr(aTmp,'bestNFeatures')),'in_canalysis'] <- TRUE
      # PREPARE CDATA
      cdata <- set_cdata(m,settings=settings, prefix=paste2('With_',i,'_Features'))
      # PREPARE CRESULT
      x[[as.character(i)]] <- set_cresult(cdata=cdata,
         fun_stats=list(), fun_plot=list(), fun_pattern=NULL, fun_bic_pattern=NULL,
         nbr_top_models=5, cfun=fun_mbc_em,
         cfun_settings=list(modelName=clusterModel, G=clusterG, rseed=clusterRseed))
      x[[as.character(i)]] <- fun_cmodel(x[[as.character(i)]])
   }
   # V, SCALE UP
   combi <- t(combn(names(x),2))
   compList <- list()
   v<- matrix(NA,length(x),length(x),dimnames=list(names(x),names(x)))
   for(i in 1:nrow(combi)){
       idx <- paste(combi[i,],collapse=",")
       compList[[idx]] <- compare_transform(x[[combi[i,1]]],x[[combi[i,2]]])
       v[combi[i,2],combi[i,1]] <- v[combi[i,1],combi[i,2]] <- compList[[idx]]$summary[1]
    }
    print(v)
    # GRAPHIC OUTPUT
   fileName <- paste2(today(),'analysisVarFeat')
   if(postscript){
      postscript(file=paste2(fileName,'.ps'))
      image(1L:nrow(v), 1L:ncol(v), v, xlim = 0.5 + c(0, nrow(v)),ylim = 0.5 +
         c(0, ncol(v)),    axes = FALSE, xlab = "", ylab = "",
         col=brewer.pal(9,"Greys"), main='scaling up
         (V)',zlim=range(v,na.rm=TRUE))
      fv_idx <- 1:ncol(v)
      if(ncol(v) > 20) 
         fv_idx <- as.integer(quantile(fv_idx))
         axis(2, at=(1L:ncol(v))[fv_idx], labels = colnames(v)[fv_idx], las = 2, line = -0.5, tick = 0)
         axis(1, 1L:nrow(v), labels = row.names(v),las = 1, line = -0.5, tick = 0)
         graphics.off()
   }
   # SAVE
   cat("\nSave 'cresult' into ", fileName)
   save(list = "x", file = paste2(fileName,'.RData'))
   return(x)
}

