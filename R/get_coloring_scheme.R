`get_coloring_scheme` <-
function(data, canalysis_variables = parkinson_canalysis_variables){
   data <- data.matrix(data)
   # MAKE A DENDROGRAM FROM THESE AVG PATTERN FOR THE DENDROGRAM ORDERING
   dendro_patterns      <- hclust(dist(data))
   # SELECT DIVERGING COLORS FROM RColorBrewer
   cluster_colors       <- brewer.pal(nrow(data),"Set1")
   # RETURN A DATA STRUCTURE WITH CLUSTERS RIGHTLY ORDERED AND THEIR ASSOCIATED COLOR
   return(list(cluster_order=dendro_patterns$order,cluster_color=cluster_colors))
}

