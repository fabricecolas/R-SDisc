\name{model_knn}
\alias{model_knn}
\title{ Returns a k-NN prediction (model) given a training and test set}
\description{ Returns a k-NN prediction (model) given a training and test set. 

Used as a generic interface to the different machine learning algorithms.}
\usage{ model_knn(formula, trainset, testset, k = 1, l = 0, prob = FALSE,
use.all = TRUE) }
\arguments{
  \item{formula}{ unused,}
  \item{trainset}{ the set of observations used for training, }
  \item{testset}{ the set of observations used for testing, }
  \item{k}{ the number \code{k} of nearest neighbours, }
  \item{l}{ the minimum number of votes for definite decision, }
  \item{prob}{ if true, then the proportion of the votes for the winning class
  are returned as attribute 'prob',}
  \item{use.all}{ controls handling of ties. If true, all distances equal to
  the 'k'th largest are included. If false, a random selection of distances
  equal to the 'k'th is chosen to use exactly 'k' neighbours.}
}
\seealso{
\code{\link{model_knn}},
\code{\link{model_naive_bayes}},
\code{\link{model_svm}},
\code{\link{predict_knn}},
\code{\link{predict_naive_bayes}},
\code{\link{predict_svm}},
\code{\link[e1071:e1071-package]{knn}},
\code{\link[e1071:e1071-package]{naiveBayes}}, 
\code{\link[e1071:e1071-package]{svm}}}
\author{ Fabrice Colas }
