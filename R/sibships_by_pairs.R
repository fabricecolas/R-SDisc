`sibships_by_pairs` <-
function(data){
	# REMOVE INDIVIDUALS IN SIBSHIPS >= 3
	individual_list <- dimnames(data[data[,"family"]  %in% names(table(data[,"family"])[table(data[,"family"])>=2]) ,])[[1]]
	# FOR EACH FAMILY (ROW), LIST THE ASCERTAINED SIBLINGS (COLUMNS: 1,2,3,4,5)
	sibling_list		<- reshape(data[individual_list,c("UNIEK","family","member")],
							timevar="member", idvar="family", direction="wide")
	siblingpair_list	<- c()
	# FOR EACH FAMILY, CONSIDER ONLY THE TWO FIRST SIBLINGS
	for(s in 1:nrow(sibling_list))
		siblingpair_list<- rbind(siblingpair_list,t(sibling_list[s,!is.na(sibling_list[s,])][2:3]))
	return(list(data=data[siblingpair_list,],model=NULL,type="sibpairs"))
}

