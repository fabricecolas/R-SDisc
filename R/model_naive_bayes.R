`model_naive_bayes` <-
function(formula, trainset, testset, laplace = 1){
   return(naiveBayes(class ~ . , data=trainset, laplace = 1))
}

