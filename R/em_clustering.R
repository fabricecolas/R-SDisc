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

