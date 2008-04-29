\name{graphic_histograms}
\alias{graphic_histograms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphic_histograms(cdata, fileName = "2007-11-24_IMG_distributions.ps", xlim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdata}{ ~~Describe \code{cdata} here~~ }
  \item{fileName}{ ~~Describe \code{fileName} here~~ }
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
function(cdata,fileName="2007-11-24_IMG_distributions.ps",xlim=NULL){
        local_data <- cdata[["data"]][,cdata[["canalysis_variables"]]]
        local_data <- as.data.frame(local_data)
	par('mai'=c(0.5,0.5,0.5,1),'mfrow'=c(4,5))
        #
        hist_list  <- list()
        hist_stats <- list()
        hist_counts <- c()
        models <- cdata[["model"]]
        for(i in colnames(local_data)){
           hist_list[[i]]  <- hist(local_data[,i],plot=FALSE)
           local_stat_txt <- ""
           hist_stats[[i]] <- list(txt="",xlim=list(),ylim=list())
           for(j in names(models)){
               local_model <- models[[j]][["model"]][[i]]
               if(models[[j]][["type"]] == "lm")
                  local_stat_txt <- paste2(local_stat_txt,local_model[["formula"]],"\n")
               else
                  local_stat_txt <- paste2(local_stat_txt,sprintf("\%1.2f",local_model[[1]]),", ")
           }
           local_stat_txt <- paste2(local_stat_txt,length(table(local_data[,i]))," distinct values")
           hist_stats[[i]][["txt"]] <- local_stat_txt
           hist_stats[[i]][["xlim"]] <- c(min(local_data[,i],na.rm=TRUE),max(local_data[,i],na.rm=TRUE))
#           hist_stats[[i]][["xlim"]] <- c(quantile(data.matrix(local_data[,i]),probs=c(0.05),na.rm=TRUE),
#                quantile(data.matrix(local_data[,i]),probs=c(0.95),na.rm=TRUE))
           hist_stats[[i]][["ylim"]] <- range(hist_list[[i]][["counts"]])
        }
        #
	for(i in 1:ncol(local_data)){
            par('lab'=c(3,3,3),  # number of ticks
                'mar'=c(2.5,2.5,1.5,1) # margins bottom,left,top,right
                )
            plot(hist_list[[i]],
                main=dimnames(local_data)[[2]][i],
                xlab="",ylab="",xlim=hist_stats[[i]][["xlim"]],ylim=hist_stats[[i]][["ylim"]]
                #,ylog=TRUE
                ) #!!!!
            text(x=mean(hist_stats[[i]][["xlim"]]),y=hist_stats[[i]][["ylim"]][2]*0.85,hist_stats[[i]][["txt"]])
            }
		
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
