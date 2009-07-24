\name{stratified_traintest_split}
\alias{stratified_traintest_split}
\title{Returns a list of index vectors with training/test split}
\description{ Stratified training test sampling to train machine learning
classification algorithms. The sampling strategy preserves in each single
training/test set the cluster membership proportion. Therefore, we discard from
the observed variability, the uncertainty associated to the cluster memberships
proportions.}
\usage{ stratified_traintest_split(data, K = 10, perc = 0.7, rseed = 6013) }
\arguments{
  \item{data}{ A data matrix, }
  \item{K}{ Number of random draws, }
  \item{perc}{ The percentage of observation in the training set (\code{perc})
      and in the test set (\code{perc}), }
  \item{rseed}{ For reproducibility of the result, we set he random seed.}
}
\author{ Fabrice Colas }
