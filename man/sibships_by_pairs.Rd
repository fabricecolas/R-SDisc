\name{sibships_by_pairs}
\alias{sibships_by_pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sibships_by_pairs(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ ~~Describe \code{data} here~~ }
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
function(data){
	# REMOVE INDIVIDUALS IN SIBSHIPS >= 3
	individual_list <- dimnames(data[data[,"family"]  \%in\% names(table(data[,"family"])[table(data[,"family"])>=2]) ,])[[1]]
	# FOR EACH FAMILY (ROW), LIST THE ASCERTAINED SIBLINGS (COLUMNS: 1,2,3,4,5)
	sibling_list		<- reshape(data[individual_list,c("UNIEK","family","member")],
							timevar="member", idvar="family", direction="wide")
	siblingpair_list	<- c()
	# FOR EACH FAMILY, CONSIDER ONLY THE TWO FIRST SIBLINGS
	for(s in 1:nrow(sibling_list))
		siblingpair_list<- rbind(siblingpair_list,t(sibling_list[s,!is.na(sibling_list[s,])][2:3]))
	return(list(data=data[siblingpair_list,],model=NULL,type="sibpairs"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
