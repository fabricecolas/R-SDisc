\name{analysis}
\alias{analysis}
\title{Sample workflow-method for subtype discovery}
\description{ Workflow example for subtype discovery by cluster analysis. This
workflow generates the graphics that illustrate the selected data
transformation.  Then, it proceeds to a cluster analysis. Further, we
characterize the cluster results by table-reports. And we finally cross-compare
the top most likely results.  }
\usage{ analysis(x, device="PS", img=FALSE, html=TRUE) }
\arguments{ 
\item{x}{Of type \code{cresult} which can be constructed from the
method \code{set_cresult}, } 
\item{device}{ A character string telling the type of device to output the
plots (in \code{"PS", "PNG"}). }
\item{img}{A logical telling whether the HTML report should also include images,  } 
\item{html}{A logical telling whether an HTML report should be generated.  } }
\author{ Fabrice Colas }
\keyword{ cluster }
\keyword{ models }
\keyword{ design }
\keyword{ manip }
\keyword{ hplot }
