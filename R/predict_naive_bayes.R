`predict_naive_bayes` <-
function(model,testset){
   return(map(predict(model,testset[,-match("class",names(testset))],type="raw")))
}

