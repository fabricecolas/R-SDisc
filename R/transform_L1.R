`transform_L1` <-
function(data){
	data <- data/matrix(apply(data,2,sum,na.rm=TRUE),nrow(data),ncol(data),byrow=TRUE)
	# TODO: UPDATE TYPE FIELD	
	return(list(data=data,model=NULL,type="L1"))
}

