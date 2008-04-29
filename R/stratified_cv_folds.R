`stratified_cv_folds` <-
function(n, K, G) 
{
   set.seed(6013)
   # WE DO A ROTATION ON THE ROWS OF THE FINAL FOLD MATRIX FOR NOT HAVING
   # ALWAYS THE SAME FOLD-STRATAS WITH ONE-LESS COUNT, I.E. THE ONE FOLDS AT
   # THE END OF THE MATRIX, 
   if(G == K) 
      G <- 1
   rowPermutation <- c(((G %% K)+1):K,1:(G %% K))
   # 
   size <- n/K
   cv <- matrix(0, K, ceiling(size))
   if (size < 5 & size != 1){
      print("The number of folds is too large for the data size")
      return(NULL)
   }
   if (size == 1) 
      return(matrix(1:n, nrow = n))
   size.int <- floor(size)
   size.vector <- rep(size.int, K)
   if (size.int != size) 
      size.vector[1:((size - size.int) * K)] <- size.vector[1:((size - size.int) * K)] + 1
   group.index <- c()
   for (j in 1:K) 
      group.index <- c(group.index, rep(j, size.vector[j]))
   group.index <- group.index[sample(n, n, replace = FALSE)]
   for (j in 1:K) {
      whichj <- which(group.index == j)
      cv[j, 1:length(whichj)] <- whichj
   }
   return(cv[rowPermutation,])
}

