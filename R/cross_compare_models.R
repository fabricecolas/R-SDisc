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

