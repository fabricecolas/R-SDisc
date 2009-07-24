\name{model_svm}
\alias{model_svm}
\title{ Returns a SVM model given a training set}
\description{ Returns a SVM model given a training set. 

Used as a generic interface to the different machine learning algorithms.}
\usage{ model_svm(formula, trainset, testset, kernel = "linear", type =
"C-classification") }
\arguments{
  \item{formula}{ the formula for the training of the classifier, }
  \item{trainset}{ the data considered as trainset, }
  \item{testset}{ unused,}
  \item{kernel}{ the kernel used in training and predicting. You might
  consider changing some of the following parameters, depending on the kernel
  type.  linear, polynomial, radial basis, sigmoid.}
  \item{type}{ \code{svm} can be used as a classification machine, as a
  regresson machine, or for novelty detection. Depending of whether 'y' is a
  factor or not, the default setting for 'type' is 'C-classification' or
  'eps-regression', respectively, but may be overwritten by setting an explicit
  value.  Valid options are: 'C-classification' 'nu-classification'
  'one-classification' (for novelty detection) 'eps-regression'
  'nu-regression'.  }
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
