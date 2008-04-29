`statistics_logodds` <-
function(data,class,fun_midthreshold=median){
        #
        class           <- map(class)
	local_data      <- cbind(data[["transformed"]],class=class)
	s               <- as.data.frame(matrix(0,length(unique(local_data$class)),length(data[["sumscore_groups"]])))
	# FOR EACH HIERARCHICAL SUBSET OF OUTCOMES, MAKE SUM SCORES
	for(gr in names(data[["sumscore_groups"]])){
		mat.l   <- cbind(SScore=apply(data.matrix(local_data[,data[["sumscore_groups"]][[gr]]]),1,sum),class=local_data$class)
		# FOR EACH PATTERN, DERIVE ITS STATISTIC
		for(g in sort(unique(local_data$class))){
			med_SScore <- fun_midthreshold(mat.l[,"SScore"])
			m11 <- nrow(mat.l[mat.l[,"class"]==g & mat.l[,"SScore"] >= med_SScore,])
			if(is.null(m11))
				m11 <- 0
			m12 <- nrow(mat.l[mat.l[,"class"]==g & mat.l[,"SScore"] < med_SScore,])
			if(is.null(m12))
				m12 <- 0			
			m21 <- nrow(mat.l[mat.l[,"class"]!=g & mat.l[,"SScore"] >= med_SScore,])
			if(is.null(m21))
				m21 <- 0			
			m22 <- nrow(mat.l[mat.l[,"class"]!=g & mat.l[,"SScore"] < med_SScore,])
			if(is.null(m22))
				m22 <- 0			
			mat <- matrix(c(m11,m12,m21,m22),2,2,byrow=TRUE)
			# LOG OF THE ODD-RATIO ALSO NAMED CROSS PRODUCT
			s[g,match(gr,names(data[["sumscore_groups"]]))] <- log(mat[1,1] * mat[2,2] / (mat[1,2] * mat[2,1]))
			dimnames(s)[[2]][match(gr,names(data[["sumscore_groups"]]))] <- paste(gr,sprintf("%.1f",med_SScore),sep="_",collapse="")
		}
	}
	return(list(out=s))
}

