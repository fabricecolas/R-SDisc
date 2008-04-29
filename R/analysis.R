`analysis` <-
function(data,
	canalysis_variables = list(),     # OUTCOMES ON WHICH THE TRANSFORMATION AND CLUSTERING SHOULD BE DONE
	sumscore_groups = list(),       # LOGODDS AND EVENTUALLY OTHER STATISTICS USE THIS GROUPING TO DO SUM SCORES
        stats_fun = list(),             # A LIST OF FUNCTION TO EVALUATE THE CLUSTERING RESULT
        G_set=1:9,
        modelNames         = c("EII","VII","EEI","EVI","VEI","VVI","EEE","EEV","VEV","VVV"),
        formula_matrix = NULL,
	fun_transform=list(identity=function(x){return(data=I(x),model=NULL)}),filePrefix=NULL,...){
	#
	# EVENTUALLY INITIALIZE THE FIRST PART OF THE FILE-NAMING 
	# AND THE LIST OF OUTCOMES ON WHICH THE ANALYSIS WILL BE PERFORMED
	#
	if(!is.null(filePrefix)){
            filePrefix <- paste2(today(),"_",filePrefix)
            postscript(file=paste2(filePrefix,"Cluster_Analysis.ps"))
            }
	if(is.null(canalysis_variables))
            canalysis_variables <- dimnames(data)[[2]]
	#
	# 0. DATA SELECTION
	#		CONFINE TO COMPLETE CASES AND PROCEED TO THE ANALYSIS
	#
	data <- list(   orig    = data,
                        cc      = data[row.names(na.omit(data[,canalysis_variables])),],
                        canalysis_variables  = sort(canalysis_variables),
                        sumscore_groups    = sumscore_groups,
                        formula_matrix     = formula_matrix,
                        mask_canalysis     = colnames(formula_matrix) %in% canalysis_variables,
                        mask_evaluation    = colnames(formula_matrix) %in% unlist(sumscore_groups))
	#
	# 0. DATA SELECTION
	#		PROCEED TO AN ANALYSIS ON THE OUTCOMES 
	#	. CORRELATION
	#	. BOXPLOTS
	#	. HISTOGRAMS
	# 
        orig_data <- list(data=data[["cc"]],model=NULL,canalysis_variables=canalysis_variables)
        variables_analysis(orig_data,filePrefix=paste2(filePrefix,"Original_Data_CC"))
	#
	# 1. AND 2. MULTIPLE CLUSTERING ANALYSIS AND MODEL SELECTION (BIC)
	#	. MODEL BASED CLUSTERING FROM 1-9 GROUPS WITH DIFFERENT COVARIANCE MATRICES
	#	. BIC TABLE FOR MODEL AND NUMBER OF GROUP SELECTION
	#	. DATA TRANSFORMATION
	# 
        transform_out <- transform_variables(data,fun_transform=fun_transform)
	data[["transformed"]] <- transform_out[["data"]]
	data[["fun_transform"]] <- fun_transform
	#
	# ALSO ANALYSE OUTCOMES FOR TRANSFORMED DATA
	#
	# variables_analysis(data[["transformed"]][,canalysis_variables],filePrefix=paste2(filePrefix,"Transformed_Data"))
        transform_out[["canalysis_variables"]] <- canalysis_variables
        variables_analysis(transform_out,filePrefix=paste2(filePrefix,"Transformed_Data"))
	#
	canalysis               <- cluster_analysis(data,stats_fun = stats_fun, filePrefix=filePrefix,G_set=G_set,modelNames=modelNames)
	data[["cc"]]            <- cbind(data[["cc"]],class=map(canalysis$z))
	data[["transformed"]]   <- cbind(data[["transformed"]],class=map(canalysis$z))
	#
	# CHARACTERIZE EACH CLUSTER BY ITS AVG PATTERN 
        #
	#
	# 3. ASSESS THE CLUSTERING METHOD CONSISTENCY
	#	DO MODEL COMPARISON ROW- AND COLUMN-WISE IN THE 5% BOUNDING BOX 
	# 
        data[["cluster_analysis"]] <- canalysis
	compare_canalysis(data,filePrefix=filePrefix,bbox_threshold=2)
        save_models(data,filePrefix=filePrefix)
	#
	# 4. CHARACTERIZE BY VISUALIZATION THE DIFFERENT CLUSTERS
	#	. PARALLEL COORDINATES (USER-DEFINED)
	#	. HEATMAPS (DEFAULT)
	#	. DENDROGRAMS ILLUSTRATING AVG-PATTERNS AND OUTCOME'S SIMILARITY
	# 
	if(!is.null(filePrefix))
                graphics.off()
	return(data)
}

