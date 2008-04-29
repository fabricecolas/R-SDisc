`compare_canalysis` <-
function(analysis1,analysis2=NULL,query=NULL,filePrefix=NULL,bbox_threshold = 5){
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

