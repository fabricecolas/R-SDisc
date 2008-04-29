`transform_ALL` <-
function(data,model=NULL,type="min",effect="duration"){
	if(is.null(model)){
		run.vector <- dimnames(data)[[2]]
		model <- list()
		for(i in run.vector){
			if(type == "max"){
				model[[i]] <- max(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "sigma"){
				model[[i]] <- sd(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "absmax"){
				model[[i]] <- max(abs(data[,i],na.rm=TRUE),na.rm=TRUE)
				data[,i] <- data[,i]/model[[i]]
				}
			if(type == "min"){
				model[[i]] <- min(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			if(type == "avg"){
				model[[i]] <- mean(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			if(type == "median"){
				model[[i]] <- median(data[,i],na.rm=TRUE)
				data[,i] <- data[,i]-model[[i]]
				}
			}
		return(list(data=data,type=type,model=model))
	}
	else{
		for(o in names(model)){
			if(type == "max" || type == "absmax" || type == "sigma")
				data[,o] <- data[,o]/model[[o]]
			if(type == "min" || type == "median" || type == "avg")
				data[,o] <- data[,o]-model[[o]]				
		}
		return(list(data=data,type=type,model=model))
	}
}

