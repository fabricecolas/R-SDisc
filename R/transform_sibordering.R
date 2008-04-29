`transform_sibordering` <-
function(data){
   for(f in unique(data$family)){
      fmembers <- sort(data[data$family == f,"member"])
      if(!is.na(fmembers[1]) && !is.na(fmembers[2])){
         data[data$family == f & data$member == fmembers[1],"member"] <- 1
         data[data$family == f & data$member == fmembers[2],"member"] <- 2 
      }
   }
   return(data[data$member <= 2,])
}

