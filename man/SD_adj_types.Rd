\name{SD_adj_types}
\alias{SD_adj_types}
\docType{data}
\title{list specifying the left/right hand side for the adjustment formulas}
\description{
This list contains five elements: sqrt, log, none, square, exp. Each element
specifies the left and the right hand side part for linear models. We may use
this structure when testing and comparing different types of adjustment. 
}
\usage{data(SD_adj_types)}
\format{
  The format is:
List of 5
 $ sqrt  :List of 2
  ..$ ls: chr "sqrt("
  ..$ rs: chr ")"
 $ log   :List of 2
  ..$ ls: chr "log("
  ..$ rs: chr ")"
 $ none  :List of 2
  ..$ ls: chr ""
  ..$ rs: chr ""
 $ square:List of 2
  ..$ ls: chr "("
  ..$ rs: chr ")^2"
 $ exp   :List of 2
  ..$ ls: chr "exp("
  ..$ rs: chr ")"
}
\details{
}
\source{
  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
\emph{Stability of Clusters for Different Time Adjustments in Complex Disease Research}.
under review for IEEE EMBC 2008.

}
\examples{
data(SD_adj_types)
## maybe str(SD_adj_types) ; plot(SD_adj_types) ...
}
\keyword{datasets}
