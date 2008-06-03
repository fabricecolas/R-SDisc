\name{fun_mbc_em}
\alias{fun_mbc_em}
\title{ Model based clustering function staring with an E-step (~ M-step?)}
\description{ This function proceeds to a model based clustering on the
provided data (a data matrix), given the model (\code{modelName}) and the
number of clusters (\code{G}). The initialization is particular in that it
draws at random a cluster membership probability matrix (\code{z}) of dimension
(N x G), with N the number of observations. From that vector it estimates the
corresponding model by an M-step and then starts EM given that model.  }
\usage{ fun_mbc_em(data, G = 3, modelName = "EII", rseed = 6013) }
\arguments{
  \item{data}{ A regular data matrix (without missing values),}
  \item{G}{ A number of clusters,}
  \item{modelName}{ A three character element referring to the type of model
  that we search for,}
  \item{rseed}{ By default, a random seed integer which we fix in order to
  initialize EM in the same way.}
}
\value{ It returns a model based clustering result.  }
\author{ Fabrice Colas }
\seealso{\link[mclust:mclust-package]{mstep}, \link[mclust:mclust-package]{em},
\link[mclust:mclust-package]{bic}, \link[mclust:mclust-package]{map}}
\keyword{}
