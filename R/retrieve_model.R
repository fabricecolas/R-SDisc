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

