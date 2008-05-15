\name{set_canalysis_result}
\alias{set_canalysis_result}
\title{ Prepare a cluster analysis result object }
\description{
   This constructor simply returns a data structure embedding the cluster result
   together with the data (\code{canalysis_data}) and the global configuration
   (\code{canalysis_config}). Its has for main element the cluster result.  
}
\usage{ set_canalysis_result(data = NULL, config = NULL, result = list()) }
\arguments{
  \item{data}{ Of type \code{canalysis_data},}
  \item{config}{ Of type \code{canalysis_config},}
  \item{result}{ A list containing the clustering algorithm output. }
}
\value{
  A data structure with four elements:
  \item{result}{ The clustering algorithm output,}
  \item{data}{ Of type \code{canalysis_data},}
  \item{config}{ Of type \code{canalysis_config},}
  \item{class}{ \code{canalysis_result}.}
}
\author{ Fabrice Colas }
\seealso{ \code{\link{set_calgo}}, \code{\link{set_canalysis_config}},
\code{\link{set_canalysis_data}}, \code{\link{set_canalysis_result}},
\code{\link{set_plot_data}}, \code{\link{set_tdata}} }
%\keyword{ ~kwd1 }
