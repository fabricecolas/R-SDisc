`cluster_analysis` <-
function(data,
        filePrefix=NULL,
        G_set=1:9,
        modelNames=c("EII", "VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),
	stats_fun = list(       oddratios=function(data,class){return(statistics_logodds(data,class))}
#				,gen_naive_bayes=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_naive_bayes,predict=predict_naive_bayes),K=5))}
#				,gen_knn=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_knn,predict=predict_knn),K=5))}
#				,gen_svm=function(data,class){return(evaluate_generalization(data,class,fun_classifier=list(model=model_svm,predict=predict_svm),K=5))}
				)){
	# 
	# PROCEED TO MBC ANALYSIS ON THE TRANSFORMED DATA
	# 
	clustering <- model_based_clustering(data,filePrefix=filePrefix,stats_fun=stats_fun,G_set=G_set)	
	return(clustering)
}

