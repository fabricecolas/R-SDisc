`compute_pattern` <-
function(data,fun=mean){
	data <- data.matrix(aggregate(data,fun,by=list(class=data$class),na.rm=TRUE))[,-1]
	row.names(data) <- 1:nrow(data)	
	if(!is.na(match("class", dimnames(data)[[2]])))
		data <- data[,-match("class", dimnames(data)[[2]])]
	return(data)
}

