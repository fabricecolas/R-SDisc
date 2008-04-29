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

