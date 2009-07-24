\name{generate_cdata_settings}
\alias{generate_cdata_settings}
\title{ Given a data set, generates a setting matrix with default values}
\description{ This function generates a settings matrix with default values. It
should be edited. This matrix is used to define the transformations to prepare
a data set. The \code{cdata} structure embeds it into its
\code{attr(,"settings")} argument.}
\usage{ generate_cdata_settings(data) }
\arguments{ \item{data}{ A regular data matrix. } }
\value{ Returns a character matrix.  }
\author{ Fabrice Colas }
\examples{
## Not run:
# LOAD DATASET
data(wada2008)
# GENERATE A SAMPLE SETTING FILE
settings_tmp<- generate_cdata_settings(wada2008)
# SHOULD BE EDITED IN A REGULAR SPREADSHEET
write.csv(settings_tmp,file="settings_tmp.csv")
settings_tmp <- read.csv('settings_tmp.csv',row.names=1)
# E.G. TO COMPARE WITH THE ONE OF WADA 2008
data(wada2008_settings)
print(wada2008_settings)
## End(Not run)
}
\seealso{ \code{\link{wada2008_settings}}, \code{\link{wada2008}}, \code{\link{set_cdata}}}
