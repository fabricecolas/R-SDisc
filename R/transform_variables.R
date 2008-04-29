`transform_variables` <-
function(data,fun_transform = NULL){
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

