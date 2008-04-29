`variables_analysis` <-
function(cdata,model=NULL,canalysis_variables,filePrefix=today(),xlim=c(-3,3)){
        local_data <- cdata[["data"]]
	#
	# CORRELATION
	#	. CORRELATION MATRIX
	#	. ABSOLUTE VALUED CORRELATION MATRIX TAKEN AS A DISTANCE MATRIX
	#	TO PLOT A DENDROGRAM WHICH ILLUSTRATES OUTCOME'S RELATION
	#
	cor_mat <- cor(local_data[,cdata[["canalysis_variables"]]],use="pairwise.complete.obs")
	write.table(cor_mat,file=paste2(filePrefix,"Correlation_Matrix.csv"),dec=",",sep=";")
	hclust(as.dist(abs(cor_mat)))
	#
	# BOXPLOTS 
	#		SUMMARY DITRIBUTIONAL SHAPE OF THE OUTCOMES BY THEIR BOX PLOTS
	#	. 1ST AND 3RD QUARTILES
	#	. MEDIAN
	#	. 95% BOUNDS
	#	. + EXTREMES BEYOND THE 95%
	# 
	par(mfrow=c(1,1),'mai'=c(1.2,0.7,0.7,0.7),new=FALSE)
	boxplot(local_data,'las'=2,main=paste2(filePrefix,"Boxplot"))
	#
	# HISTOGRAMS 
	#
	graphic_histograms(cdata,model,xlim=xlim)
}

