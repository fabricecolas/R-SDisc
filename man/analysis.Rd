\name{analysis}
\alias{analysis}
\title{Sample workflow-method for subtype discovery}
\description{ Workflow example for subtype discovery by cluster analysis. This
workflow generates the graphics that illustrate the selected data
transformation.  Then, it proceeds to a cluster analysis. Further, we
characterize the cluster results by table-reports. And we finally cross-compare
the top most likely results.  }
\usage{ analysis(cresult, doPSOutput=TRUE, csv=TRUE, html = TRUE) }
\arguments{ 
\item{cresult}{Of type \code{cresult} which can be constructed from the
method \code{set_cresult}, } 
\item{doPSOutput}{A logicial telling whether doPSOutput output should be
generated in the making of the analysis,} 
\item{csv}{A logical telling whether Excel-readable files should be prepapred
in the making of the analysis.  } 
\item{html}{A logical telling whether an HTML output should be generated.  } 
}
\author{ Fabrice Colas }
\keyword{ cluster }
\keyword{ models }
\keyword{ design }
\keyword{ manip }
\keyword{ hplot }
