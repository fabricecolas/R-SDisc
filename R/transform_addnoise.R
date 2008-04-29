`transform_addnoise` <-
function(
        data,
        model=NULL,
        canalysis_variables = NULL,
        type="addnoise",
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
   return(list(data=data,type=type,model=model))
}

