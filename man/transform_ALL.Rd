\name{transform_ALL}
\alias{transform_ALL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_ALL(data, model = NULL, type = "min", effect = "duration")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{effect}{ ~~Describe \code{effect} here~~ }
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
