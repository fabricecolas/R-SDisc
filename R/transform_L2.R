`transform_L2` <-
function(data){
	data <- data/sqrt(matrix(apply(data^2,2,sum,na.rm=TRUE),nrow(data),ncol(data),byrow=TRUE))
	# TODO: UPDATE TYPE FIELD
	return(list(data=data,model=NULL,type="L2"))
}

