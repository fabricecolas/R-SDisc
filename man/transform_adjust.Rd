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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
