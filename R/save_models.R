`save_models` <-
function(data,filePrefix=NULL,threshold=5){
   m_set                <- get_model_set(data[["cluster_analysis"]][["bic_table_relative"]],bbox_threshold=threshold)
   for(l in 1:length(data[["cluster_analysis"]][["out"]])){
      local_canalysis   <- data[["cluster_analysis"]][["out"]][[l]]
      g_match           <- match(local_canalysis[["G"]],m_set[,"G"])
      m_match           <- match(local_canalysis[["modelName"]],m_set[,"modelName"])
      if( !is.na(g_match) & !is.na(m_match)){
	data_out        <- cbind(         ids=row.names(data[["transformed"]])
                                        , class=map(local_canalysis[["z"]])
                                        , local_canalysis[["z"]])
        colnames(data_out)[3:ncol(data_out)] <- 1:local_canalysis[["G"]]
        fileNameCSV     <- paste2(filePrefix,"_Class_",local_canalysis[["modelName"]],"-",local_canalysis[["G"]],".csv") 
        print(fileNameCSV)
        write.csv2(data_out,file=fileNameCSV)
      }
   }
}

