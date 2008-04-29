`get_vertice_set` <-
function(relativeBic, bbox_threshold=5){
   relativeBic <- list( bicTable=relativeBic, 
                        # TO SELECT THE OPTIMAL ONES BY RANKING, I.E. THE 5 BEST ONES 
                        bbox = apply(relativeBic <= max(sort(relativeBic)[1:bbox_threshold]), 2,which)
                        # TO SELECT ALL MODELS BELOW A GIVEN THRESHOLD
                        # bbox = apply(abs(relativeBic) < bbox_threshold,2,which)
                        )
   relativeBic$bbox <- relativeBic$bbox[which(lapply(relativeBic$bbox,length) >0)]
   # DERIVE VERTICES TO COMPARE WITHIN THE 5% BBOX OF THE BIC TABLE
   relativeBic$vertices <- matrix(0,0,2,dimnames=list(list(),list("m1_g1","m2_g2")))
   for(m1 in names(relativeBic$bbox))
      for(m2 in names(relativeBic$bbox))
         if(match(m2,names(relativeBic$bbox)) >= match(m1,names(relativeBic$bbox)))
            for(g1 in relativeBic$bbox[[m1]])
               for(g2 in relativeBic$bbox[[m2]])
                  if((m1 != m2 & g2 >= g1) | (m1 == m2 & g2 > g1))
                     relativeBic$vertices <- rbind(relativeBic$vertices,c(
                                                      paste(m1,row.names(relativeBic[["bicTable"]])[g1],collapse="",sep=","),
                                                      paste(m2,row.names(relativeBic[["bicTable"]])[g2],collapse="",sep=",")))
   return(relativeBic[["vertices"]])
}

