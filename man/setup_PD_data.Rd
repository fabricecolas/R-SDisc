\name{setup_PD_data}
\alias{setup_PD_data}
\docType{data}
\title{Data setup file for the Parkinson Disease analysis}
\description{ Given a data set, a basic settings file can be generated from
\code{generate_cdata_settings} and it should be edited then.}
\usage{data(setup_PD_data)}
\format{
  A data frame with 20 observations on the following 5 variables.
  \describe{
    \item{\code{group}}{a factor with levels \code{autonomic} \code{beck} \code{cognition} \code{duration} \code{dyskinesia} \code{eds} \code{motfluct} \code{motricity} \code{psychosis} \code{sleep}}
    \item{\code{in_canalysis}}{a logical vector}
    \item{\code{fun_transform}}{a factor with levels \code{transform_AVG
    transform_SIGMA} \code{transform_adjust('autcardio~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('autgi~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('auturin~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('beck~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('cogattention~duration')
    transform_AVG transform_SIGMA}
    \code{transform_adjust('cogexecutief~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('cogmemory~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('cogvisuo~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('dysk~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('eds~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('motfluct~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('psychosis~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('sleep~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('spesaxial2~duration')
    transform_AVG transform_SIGMA} \code{transform_adjust('spesaxial~duration')
    transform_AVG transform_SIGMA} \code{transform_adjust('spesbrady~duration')
    transform_AVG transform_SIGMA}
    \code{transform_adjust('spesfreezing~duration') transform_AVG
    transform_SIGMA} \code{transform_adjust('spesrigiditeit~duration')
    transform_AVG transform_SIGMA}
    \code{transform_adjust('spestremor~duration') transform_AVG
    transform_SIGMA}}
    \item{\code{visu_groups}}{a factor with levels \code{Autonomic and Sumscores} \code{Cognitive and Motricity} \code{Other}}
    \item{\code{visu_ycoord}}{a numeric vector}
  }
}
\examples{ data(setup_PD_data) }
\seealso{\code{\link{generate_cdata_settings}}}
\keyword{datasets}
