\name{stats_generalization}
\alias{stats_generalization}
\title{ Estimate the generalization of a classification algorithm given on a
dataset}
\description{ We evaluate the algorithms using classification accuracy which we
estimate by rwpeating ten times on a training-test split with 70\% of the
individuals in the training set and 30\% in the test set and discarding the
variability that depend on the cluster membership probability. In every
training and test set, we do this by preserving the cluster proportions of the
original sample. From the viewpoint of machine learning, this practice is
coherent because we primarily aim to compare algorithms. But from the viewpoint
of medical research, we might prefer to alter this stratified sampling strategy
to include the variability inherent to the cluster proportions. }
\usage{ stats_generalization(data, class, fun_classifier = list(model =
model_naive_bayes, predict = predict_naive_bayes), K = 7, fun_eval =
stratified_traintest_split, title = NULL) }
\arguments{
  \item{data}{ a data set to classify of type (\code{canalysis_data}) which
  stores especially the transformation estimates, }
  \item{class}{ a labelling given a cluster model,}
  \item{fun_classifier}{ a named list of two elements defining a function to
  \code{model} and to \code{predict},}
  \item{K}{ the number of repeats, }
  \item{fun_eval}{ by default, \code{stratified_traintest_split} where cluster
  priors remain the same in all folds,}
  \item{title}{ the title to use when outputting table of results. }
}
\author{ Fabrice Colas }
\keyword{}
