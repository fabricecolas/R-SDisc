\name{variables_analysis}
\alias{variables_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
variables_analysis(cdata, model = NULL, canalysis_variables, filePrefix = today(), xlim = c(-3, 3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdata}{ ~~Describe \code{cdata} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{canalysis_variables}{ ~~Describe \code{canalysis_variables} here~~ }
  \item{filePrefix}{ ~~Describe \code{filePrefix} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
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
	#	. 95\% BOUNDS
	#	. + EXTREMES BEYOND THE 95\%
	# 
	par(mfrow=c(1,1),'mai'=c(1.2,0.7,0.7,0.7),new=FALSE)
	boxplot(local_data,'las'=2,main=paste2(filePrefix,"Boxplot"))
	#
	# HISTOGRAMS 
	#
	graphic_histograms(cdata,model,xlim=xlim)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
