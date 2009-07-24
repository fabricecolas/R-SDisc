\name{model_naive_bayes}
\alias{model_naive_bayes}
\title{ Returns a naive Bayes model given a training set}
\description{ Returns a naive Bayes model given a training set. 

Used as a generic interface to the different machine learning algorithms.}
\usage{ model_naive_bayes(formula, trainset, testset, laplace = 1) }
\arguments{
  \item{formula}{ the formula for the training of the classifier, }
  \item{trainset}{ the data considered as trainset, }
  \item{testset}{ unused,}
  \item{laplace}{ positive double controlling Laplace smoothing}
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
