\name{graphic_characterization}
\alias{graphic_characterization}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphic_characterization(data, class_count = class_count, graphics = graphic_params, titlePrefix = NULL, canalysis_variables = canalysis_variables)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{class_count}{ ~~Describe \code{class_count} here~~ }
  \item{graphics}{ ~~Describe \code{graphics} here~~ }
  \item{titlePrefix}{ ~~Describe \code{titlePrefix} here~~ }
  \item{canalysis_variables}{ ~~Describe \code{canalysis_variables} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(data, class_count=class_count,
	graphics=graphic_params,
	titlePrefix=NULL,
        canalysis_variables = canalysis_variables){
        #
        local_data <- data 
	# LOAD TRANSFORMATION MATRIX (X,Y,Z)
	visu_mat <- read.csv(file_longitudinal_visu,header=TRUE,row.names=1,dec=",")
	# visu_mat <- visu_mat[complete.cases(visu_mat[dimnames(data[[1]])[[2]],]),]
	vignettes.title	<- levels(visu_mat$Z)[levels(visu_mat$Z) != 0]
        # WATCH-OUT! THE PATTERN [[1]] IS EXPECTED TO BE THE ONE FROM WHICH THE
        # ORDERING IS DERIVED THUS, THESE SHOULD REFER TO THE MEAN OR THE
        # MEDIAN PATTERNS
        color_sel       <- get_coloring_scheme(data[[1]],canalysis_variables = canalysis_variables)
	par(mfrow=c(3,3), 'las'=1,'mai'=c(0.4,0.55,0.18,0.05))
	#
	# LOOP OVER EACH VIGNETTE
	#
	for(v in vignettes.title){
		# TAKE THE SUBSET OF PHENOTYPES TO BE VISUALIZED ON VIGNETTE Z: TRANSFORMATION MATRIX
		T_s					<- visu_mat[visu_mat$Z == v,]
		plot(x=0        ,new=TRUE
                                ,ann=FALSE
                                ,pch=18
                                ,xlim=c(graphics[["xlim"]]$min,graphics[["xlim"]]$max)
                                ,ylim=c(graphics[["ylim"]]$min,graphics[["ylim"]]$max)
                                ,col="white",axes=FALSE)
                title(main=paste2(titlePrefix,v))
		# LOOP OVER EACH GROUP PATTERN THAT IS TO VISUALIZE
		for(g in color_sel[["cluster_order"]]){
			# LOOP OVER THE DIFFERENT STATISTICS FOR EACH PATTERN
			for(s in names(data)){
                                local_pattern           <- data[[s]][color_sel[["cluster_order"]],]
				C_i_s                   <- cbind(as.numeric(local_pattern[g,dimnames(T_s)[[1]]]),0)
				D_i_s			<- C_i_s+T_s[,1:2]
				D_i_s			<- D_i_s[sort.list(D_i_s$Y),]
                                gap                     <- D_i_s[2,"Y"] - D_i_s[1,"Y"]
				for(l in 1:(nrow(D_i_s)-1)){
                                        is_white_gap <- (D_i_s[l+1,"Y"] - D_i_s[l,"Y"] > gap)
					if((s == "median" || s == "avg") & !is_white_gap )
						arrows(D_i_s[l,1],D_i_s[l,2],D_i_s[l+1,1],D_i_s[l+1,2],col=color_sel[["cluster_color"]][g],length=0,lwd=3)
                                        if((s == "lowquant" || s == "upquant") & !is_white_gap)
						arrows(D_i_s[l,1],D_i_s[l,2],D_i_s[l+1,1],D_i_s[l+1,2],col=color_sel[["cluster_color"]][g],length=0,lty="dashed",lwd=1)
                                        # else, as (is_white_gap == TRUE) then do not draw any arrow...
				}
				if(v == graphics$legend_on){
					# PLOT THE LEGEND ON ONE OF THE VIGNETTES			
					# DO THE LEGEND WITHIN [3.5;6.5] AND START FROM THE TOP, I.E. 6.5
					#y1 <- y0 <- 0.9*graphics[["ylim"]]$max-(g-1)*3/nrow(data[[1]])
					y1 <- y0 <- graphics[["legend_ycoord"]] + 5*g/nrow(data[[1]])
					x0 <- graphics[["legend_xcoord"]]
					x1 <- 1.2*graphics[["legend_xcoord"]]
					arrows(x0,y0,x1,y1,col=color_sel[["cluster_color"]][g],length=0,lwd=3)
                                        local_g <- color_sel[["cluster_order"]][g]
					text(x0+0.15,y1,labels=paste2(local_g," (",class_count[as.numeric(local_g)],")"),pos=4)
				}						
			}
			# DIMNAMES
			axis(2,at=D_i_s[,2],labels=dimnames(D_i_s)[[1]],las=2,tick=FALSE)
			axis(1,at=seq(from = graphics[["xlim"]]$min , to = graphics[["xlim"]]$max , by = ((graphics[["xlim"]]$max - graphics[["xlim"]]$min)/4)))
		}
	}
	# DEFINE COLOR GRADIENT
	gradient.colors <- graphics[["color_gradient"]] 
	# DETERMINES ROW AND COLUMN ORDERINGS OF THE AVG PATTERN
        data <- data[[1]][,canalysis_variables]
	rowD <- hclust(dist(data))
	colD <- hclust(dist(t(data)))
	# d <- data[rowD$order,colD$order]
        d <- data[rowD$order, rev(sort(colnames(data)))]
	# PRODUCE THE 3-PLOTS: IMAGE/HEATMAP, ROW AND COLUMN DENDROGRAMS
	image(1L:nrow(d), 1L:ncol(d), d, xlim = 0.5 + c(0, nrow(d)), ylim = 0.5 + c(0, ncol(d)), axes = FALSE, xlab = "", ylab = "",col=gradient.colors,main="Average pattern visualization")
	axis(2, 1L:ncol(d), labels = dimnames(d)[[2]], las = 2, line = -0.5, tick = 0, cex.axis = graphics[["cex"]])
	axis(1, 1L:nrow(d), labels = dimnames(d)[[1]], las = 1, line = -0.5, tick = 0, cex.axis = graphics[["cex"]])
	plot(rowD, axes = FALSE, yaxs = "i", main="Patterns arranged by similarity",ylab=NULL,cex=graphics[["cex"]])
	plot(colD, axes = FALSE, yaxs = "i", main="Outcomes arranged by similarity",ylab=NULL,cex=graphics[["cex"]])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
