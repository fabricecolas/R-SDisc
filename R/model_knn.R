`model_knn` <-
function(formula, trainset, testset, k = 1, l = 0, prob = FALSE, use.all = TRUE){
   return(knn(trainset[,-match("class",names(trainset))], testset[,-match("class",names(testset))] ,trainset[,"class"], k = k, l = l, prob = prob, use.all = use.all))
}

