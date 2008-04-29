\name{graphic_profile}
\alias{graphic_profile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphic_profile(data, psoutput = NULL, plot_sample = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{psoutput}{ ~~Describe \code{psoutput} here~~ }
  \item{plot_sample}{ ~~Describe \code{plot_sample} here~~ }
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
function(data,psoutput=NULL,plot_sample=20){
   if(!is.null(psoutput))
      postscript(file=psoutput,paper="a4")
   par(mfrow=c(3,3), 'las'=1,'mai'=c(0.3,0.3,0.3,0.5))
   id_set    <- sample(row.names(data))
   #
   for(o in colnames(data[,,1])){
      a <- as.data.frame(data[,o,])
      b <- reshape(a,direction="long",ids=row.names(a),times=1:dim(data)[[3]],timevar="Year",varying=list(names(a)),v.names=o)
      # INIT PLOT
      plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(min(a,na.rm=TRUE),max(a,na.rm=TRUE)),ylab=o)
      # RAW DATA
      for(i in 1:length(id_set)){
         if((length(id_set)-i) < plot_sample)
            line_col <- "blue"
         else
            line_col <- "black"
         lines(xy.coords(b[b$id == id_set[i],c("Year",o)]),type="l",col=line_col)
         }
      # STATS
      for(s in list(list(fun=function(x){return(mean(x,na.rm=TRUE))},col="red",lty="solid",lwd = 2),
                  list(fun=function(x){return(median(x,na.rm=TRUE))},col="red",lty="dashed",lwd = 1),
                  list(fun=function(x){return(quantile(x, na.rm=TRUE, probs = seq(0, 1, 0.05))["5\%"])},col="red",lty="dashed",lwd = 1),
                  list(fun=function(x){return(quantile(x, na.rm=TRUE, probs = seq(0, 1, 0.05))["95\%"])},col="red",lty="dashed",lwd = 1)
                  )){
         s[["pattern"]] <- cbind(Year=1:dim(data)[[3]],o=apply(data[,o,],2,s[["fun"]]))
         lines(xy.coords(s[["pattern"]]),type="l",col = s[["col"]],lty = s[["lty"]],lwd = s[["lwd"]])
      }
   }
   if(!is.null(psoutput))
      graphics.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
