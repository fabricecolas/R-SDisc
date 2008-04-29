\name{graphic_cluster_profiles}
\alias{graphic_cluster_profiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
graphic_cluster_profiles(data, canalysis_variables = parkinson_canalysis_variables, psoutput = NULL, fun_stats = list(f_mean = list(fun = function(x) {
    return(mean(x, na.rm = TRUE))
}, lty = "solid", lwd = 4), f_25quant = list(fun = function(x) {
    return(quantile(x, probs = 0.25, na.rm = TRUE))
}, lty = "dashed", lwd = 1), f_75quant = list(fun = function(x) {
    return(quantile(x, probs = 0.75, na.rm = TRUE))
}, lty = "dashed", lwd = 1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{canalysis_variables}{ ~~Describe \code{canalysis_variables} here~~ }
  \item{psoutput}{ ~~Describe \code{psoutput} here~~ }
  \item{fun_stats}{ ~~Describe \code{fun_stats} here~~ }
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
function(data,canalysis_variables = parkinson_canalysis_variables,
        psoutput        =NULL,
        fun_stats       =list(
#                f_median   = list(fun = function(x){return(median(x,na.rm=TRUE))},lty="solid",lwd = 4)
                 f_mean   = list(fun = function(x){return(mean(x,na.rm=TRUE))},  lty="solid",lwd = 4)
                , f_25quant = list(fun = function(x){return(quantile(x,probs=0.25,na.rm=TRUE))},  lty="dashed",lwd = 1)
                , f_75quant = list(fun = function(x){return(quantile(x,probs=0.75,na.rm=TRUE))},  lty="dashed",lwd = 1)
                ))
                {
   #
   if(!is.null(psoutput))
      postscript(file=psoutput,paper="a4")
   #
   par(mfrow=c(2,9), 'las'=1,'mai'=c(0.3,0.1,0.3,0.3),oma=c(1,1,1,2))
   avg_pattern  <- compute_pattern(as.data.frame(data[,c(canalysis_variables,"class"),1]),fun = fun_stats[[1]][[1]])
   color_sel    <- get_coloring_scheme(avg_pattern,canalysis_variables = canalysis_variables)
   ymin         <- min(quantile(data,probs=0.05,na.rm=TRUE),na.rm=TRUE)
   ymax         <- max(quantile(data,probs=0.95,na.rm=TRUE),na.rm=TRUE)
   #
   id_set    <- sample(row.names(data))
   for(o in  canalysis_variables){
      # INIT PLOT
      plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(ymin,ymax),main=o,ann=TRUE,axes=FALSE)
      # STATS
      for(cl in 1:length(color_sel[["cluster_order"]])){
         df_subset <- data[data[,"class",1] == color_sel[["cluster_order"]][cl],,]
         df_subset <- df_subset[row.names(df_subset)[!is.na(row.names(df_subset))],,]
         # RAW DATA
         a <- as.data.frame(df_subset[,o,])
         b <- reshape(a,direction="long",ids=row.names(a),times=1:dim(data)[[3]],timevar="Year",varying=list(names(a)),v.names=o)
         # RAW PROFILES
#         for(i in sample(unique(b$id))[1:5])
#            lines(xy.coords(b[b$id == i,c("Year",o)]),type="l",col=color_sel[["cluster_color"]][cl],lwd=1,lty="dotted")
         for(s in fun_stats){
            s[["pattern"]] <- cbind(Year=1:dim(data)[[3]],o=apply(df_subset[,o,],2,s[["fun"]]))
            lines(xy.coords(s[["pattern"]]),type="l",col = color_sel[["cluster_color"]][cl],lty = s[["lty"]],lwd = s[["lwd"]])
         }
         # PLOT THE LEGEND ON ONE OF THE VIGNETTES			
         # DO THE LEGEND WITHIN [3.5;6.5] AND START FROM THE TOP, I.E. 6.5
      }
      # Y-AXIS PLOT (2)
      axis(2,at=pretty(c(ymin,ymax),n=3))
   }
   #
   # LEGEND
   #
   plot(NULL,xlim=c(1,dim(data)[[3]]),ylim=c(ymin,ymax),main="legend",ann=TRUE,axes=FALSE)
   for(cl in 1:length(color_sel[["cluster_order"]])){
      yrange <- (ymax-ymin)/length(color_sel[["cluster_order"]])
      y1 <- y0 <- ymax-cl*0.5*yrange
      x0 <- 1
      x1 <- 2
      arrows(x0,y0,x1,y1,col=color_sel[["cluster_color"]][cl],length=0,lwd=3)
      text(1.15*x1,y1,labels=paste2(cl," (",nrow(data[data[,"class",1] == color_sel[["cluster_order"]][cl],,]),")"),pos=4)
   }
   #
   #
   #
   if(!is.null(psoutput))
      graphics.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
