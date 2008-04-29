`stratified_cv` <-
function(data,K=10){
   indexes <- list()
   G <- length(table(data$class))
   g <- as.numeric(names(sort(table(data$class))[1]))
   row_number <- nrow(data[data$class == g,])
   out_cv <- stratified_cv_folds(row_number,K=K,g)
   while(is.null(out_cv)){
      K <- K-1
      print(K)
      out_cv <- stratified_cv_folds(row_number,K=K,g)
   }
   for(g in 1:G)
      indexes[[g]] <- stratified_cv_folds(row_number,K,g)
   return(list(indexes=indexes,K=K))
}

