`predict_svm` <-
function(model, testset){
   predict(model,testset[,-match("class",names(testset))])
}

