\name{predict_svm}
\alias{predict_svm}
\title{ Returns a SVM prediction given a model}
\description{ Returns a SVM prediction given a \code{model}.  }
\usage{ predict_svm(model, testset) }
\arguments{
  \item{model}{ The training set referred to as \code{model},}
  \item{testset}{ The set of examples to classify referred to as the
  \code{testset}.}
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
