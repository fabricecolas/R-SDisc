`model_svm` <-
function(formula, trainset, testset, kernel="linear",type="C-classification"){
   return(svm(class ~ . , data=trainset,kernel=kernel,type=type))
}

