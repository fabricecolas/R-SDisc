\name{statistics_logodds}
\alias{statistics_logodds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
statistics_logodds(data, class, fun_midthreshold = median)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{class}{ ~~Describe \code{class} here~~ }
  \item{fun_midthreshold}{ ~~Describe \code{fun_midthreshold} here~~ }
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
			dimnames(s)[[2]][match(gr,names(data[["sumscore_groups"]]))] <- paste(gr,sprintf("\%.1f",med_SScore),sep="_",collapse="")
		}
	}
	return(list(out=s))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
