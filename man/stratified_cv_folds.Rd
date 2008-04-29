\name{stratified_cv_folds}
\alias{stratified_cv_folds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
stratified_cv_folds(n, K, G)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ ~~Describe \code{n} here~~ }
  \item{K}{ ~~Describe \code{K} here~~ }
  \item{G}{ ~~Describe \code{G} here~~ }
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
function(n, K, G) 
{
   set.seed(6013)
   # WE DO A ROTATION ON THE ROWS OF THE FINAL FOLD MATRIX FOR NOT HAVING
   # ALWAYS THE SAME FOLD-STRATAS WITH ONE-LESS COUNT, I.E. THE ONE FOLDS AT
   # THE END OF THE MATRIX, 
   if(G == K) 
      G <- 1
   rowPermutation <- c(((G \%\% K)+1):K,1:(G \%\% K))
   # 
   size <- n/K
   cv <- matrix(0, K, ceiling(size))
   if (size < 5 & size != 1){
      print("The number of folds is too large for the data size")
      return(NULL)
   }
   if (size == 1) 
      return(matrix(1:n, nrow = n))
   size.int <- floor(size)
   size.vector <- rep(size.int, K)
   if (size.int != size) 
      size.vector[1:((size - size.int) * K)] <- size.vector[1:((size - size.int) * K)] + 1
   group.index <- c()
   for (j in 1:K) 
      group.index <- c(group.index, rep(j, size.vector[j]))
   group.index <- group.index[sample(n, n, replace = FALSE)]
   for (j in 1:K) {
      whichj <- which(group.index == j)
      cv[j, 1:length(whichj)] <- whichj
   }
   return(cv[rowPermutation,])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
