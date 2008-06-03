\name{generate_canalysis_data_settings}
\alias{generate_canalysis_data_settings}
\title{ Given a data set, generates a setting matrix with default values}
\description{ This function generates a settings matrix with default values. It
should be edited. This matrix is used to define the transformations to prepare
a data set. The \code{canalysis_data} structure embeds it into its
\code{attr(,"settings")} argument.}
\usage{ generate_canalysis_data_settings(data) }
\arguments{ \item{data}{ A regular data matrix. } }
\value{ Returns a character matrix.  }
\author{ Fabrice Colas }
\seealso{ \code{\link{set_canalysis_data}}}
\keyword{}
