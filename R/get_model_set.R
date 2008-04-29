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

