\name{set_cdata}
\alias{set_cdata}
\title{ Prepare a cluster analysis data object }
\description{ This function contructs an object containing the data for a
cluster analysis. The setting table describes in particular how the data should
be processed for the cluster analysis and which are the variables used in it.
According to that information and the data set initialization procedure, the
constructor prepare the data. The main element of the resulting structure is
the transformed data set but it also data set in its original form. It also
defines the models to re-run the transformation on new data given their
modeling parameters.}
\usage{ set_cdata(data, data_o = NULL, tdata = NULL, settings,
init_fun = list(init_data_cc), prefix ="") }

\arguments{
  \item{data}{ A data set (mandatory),}
  \item{data_o}{ The original data (optional),}
  \item{tdata}{ A list of \code{tdata} (optional),}
  \item{settings}{ A matrix where we define the sequence of transformation for
  all variables (mandatory),}
  \item{init_fun}{ A \code{canalyis_data} constructor initialization function
  to run before any further processing (\code{na.omit}).}
  \item{prefix}{ A prefix which we use for all file outputs (Postscript and
  CSV). By default, we use the current date as prefix.}
}
\details{
In first place, the constructor stores a copy of the data set in \code{data_o}.

Then, it proceeds to the initialization procedure of the data set. Its default
behaviour is to discard all observations (using \code{na.omit()}) having a
missing value among the variables for the cluster analysis. 

Next, it takes the setting matrix, it reads the \code{transform} column and
from left to right, from top to bottom, it applies the transformations. To parse
the setting matrix, we first use commas and spaces, and in second the
parenthesis (for \code{transform_adjust}). 

We keep record of all models and estimated parameters so that on new data, we
can apply the same transforms. Transformation data are stored in specific
\code{tdata} structure. The \code{tdata} field of a \code{cdata} is
actually a list of list of \code{tdata} objects, whose names, at the first
list-level, refer to the variables, and at the second level, refer to the
transformation name in its sequential order of application. 
}
\value{
  \item{data}{ The transformed data set,}
  \item{data_o}{ The original data set,}
  \item{tdata}{ A list of list of \code{tdata} objects,}
  \item{settings}{ A character matrix defining the transform and the variables
  for the cluster analysis,}
 }

\author{ Fabrice Colas }
\seealso{ \code{\link{set_cdata}}, \code{\link{set_cresult}},
\code{\link{set_tdata}} }
\keyword{}
