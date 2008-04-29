`stratified_traintest_split` <-
function(data,K=10,perc=0.7){
   indexes <- list()
   G <- length(table(data$class))
   for(g in 1:G){
      strata_idx <- nrow(data[data$class == g,])
      indexes[[g]] <- matrix(0,0,strata_idx * (1-perc)) 
      set.seed(6013)
      for(k in 1:K)
            indexes[[g]] <- rbind(indexes[[g]],sample(1:strata_idx)[1:(strata_idx * (1-perc))])
      }
   return(list(indexes=indexes,K=K))
}

