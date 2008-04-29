`transform_longitudinal` <-
function(data, fun, canalysis_variables = parkinson_canalysis_variables,effect="duration"){
   data_out             <- data
   m                    <- fun(as.data.frame(data[,c(canalysis_variables,effect),1]))
   m_effect             <- 
   data_out[row.names(m[["data"]]),colnames(m[["data"]]),1]     <- data.matrix(m[["data"]])
   for(i in 2:dim(data)[[3]]){
      # WE WANT TO REMOVE THE EFFECT OF DISEASE DURATION (AKA THE TIME) FOR BASELINE,
      # BUT WE PRESERVE THE TIME EFFECT IN THE LONGITUDINAL ANALYSIS. 
      data_tmp          <- as.data.frame(cbind(data[,canalysis_variables,i],data[,effect,1]))
      colnames(data_tmp)[ncol(data_tmp)] <- effect 
      m_tmp          <- fun(data_tmp,model = m[["model"]])
      data_out[row.names(m_tmp[["data"]]),canalysis_variables,i]  <- data.matrix(m_tmp[["data"]])[,canalysis_variables]
      }
   return(data_out)
}

