\name{transform_adjust}
\alias{transform_adjust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_adjust(data, model = NULL, type = "lm", f_matrix = formula_matrix, transform = list(sqrt = list(ls = "sqrt(", rs = ")"), log = list(ls = "log(", rs = ")"), none = list(ls = "", rs = ""), square = list(ls = "(", rs = ")^2"), exp = list(ls = "exp(", rs = ")")), optimize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
  \item{model}{ ~~Describe \code{model} here~~ }
  \item{type}{ ~~Describe \code{type} here~~ }
  \item{f_matrix}{ ~~Describe \code{f_matrix} here~~ }
  \item{transform}{ ~~Describe \code{transform} here~~ }
  \item{optimize}{ ~~Describe \code{optimize} here~~ }
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
function(
        data,
        model=NULL,
        type="lm",
        f_matrix=formula_matrix,
        transform=list(
                 sqrt=  list(ls="sqrt(", rs=")")
                ,log=   list(ls="log(",  rs=")")
                ,none=  list(ls="",      rs="")
                ,square=list(ls="(",     rs=")^2")
                ,exp=   list(ls="exp(",  rs=")")
                ),
        optimize = FALSE
        )
{
   if(is.null(model)){
      # data
      data_out     <- as.data.frame(matrix(NA,nrow(data),0,dimnames=list(row.names(data),list())))
      model_out    <- list()
      # for each outcome (column of the data matrix), loop
      for(i in 1:nrow(f_matrix)){
         fact_name <- colnames(f_matrix)[i]
         # left hand side of the formula (above 0)
         lhs       <- colnames(f_matrix)[f_matrix[i,]>0] 
         # right hand side of the formula (below 0)
         rhs       <- colnames(f_matrix)[f_matrix[i,]<0]
         if(length(lhs) == 1 && length(rhs) > 0){
            list_f    <- list(paste2(paste(lhs,collapse="+"),"~"))
            # for each left hand side element
            for(f_idx in 1:length(list_f)){
               elmt <- list_f[[f_idx]]
               # we append the right hand side transform in only one way for both
               for(t_id in 1:length(transform)){
                  # in case this is the first transform, we just update list_f,
                  # otherwise we add new elements to list_f referred to as f_idx_l
                  f_idx_l <- f_idx+t_id-1
                  local_f <- ""
                  for(rhs_l in rhs){
                     local_f <- paste2(local_f,transform[[t_id]]$ls,rhs_l,transform[[t_id]]$rs)
                     if(length(rhs) > 1 && rhs_l != rhs[[length(rhs)]])
                        local_f <- paste2(local_f,"+")
                  }
                  list_f[[f_idx_l]] <- paste2(elmt,local_f)
               }
            }
            model          <- list()
            adj_rsquare    <- matrix(0,length(list_f))
            for(f_id in 1:length(list_f)){
               model[[f_id]]          <- lm(as.formula(list_f[[f_id]]),data=data)
               adj_rsquare[f_id]      <- summary(model[[f_id]])[["adj.r.squared"]]
            }
            # if there are several rsquare that are equal, take the 1st one...
            best_id                            <- which(adj_rsquare == max(adj_rsquare))[1]
            model_out[[fact_name]] <- model[[best_id]]
            model_out[[fact_name]][["rhs"]]      <- rhs
            model_out[[fact_name]][["formula"]] <- list_f[[best_id]]
            # take the corresponding residual
            r                              <- residuals(model_out[[fact_name]])
            r_names                        <- names(r)
            }# else, do not adjust data and preserve original scores
         else
            r <- NULL
         data_out                       <- cbind(data_out,matrix(NA,nrow(data_out),1))
         if(!is.null(r))
            data_out[r_names,ncol(data_out)] <- r
         else
            data_out[,ncol(data_out)]   <- data[row.names(data_out),fact_name]
         colnames(data_out)[i]          <- colnames(f_matrix)[i]
      }
      return(list(data=data_out,model=model_out,type="lm"))
   }
   else{
      for(o in names(model)){
         data[,o] <- data[,o]-predict(model[[o]],data[,c(o,unlist(model[[o]][["rhs"]]))])
         }
      return(list(data=data,type=type,model=model))
   }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
