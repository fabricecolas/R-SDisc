\name{stats_lambdasibs}
\alias{stats_lambdasibs}
\title{ Evaluate the level of familial aggregation in each cluster by the lambda sib risk ratio }
\description{ As in medical research, the clinical relevance of the clusters is
a key concern, this statistical test assesses the significance of the level of
familial aggregation in analyses involving sibling pairs. Our test relies on a
risk ratio which we refer to as the lambda sibs.

We characterize each individual as proband or sibling depending on whether this
individual was the first sibling involved in the study or not.  This test
quantifies the risk increases of the second sibling given the characteristics
of the proband. For instance, a $\lambda_{sibs} = 1$ means that the risk does
not increase and that the cluster membership of the proband does not influence
the one of his sibling. On the other hand, if $\lambda_{sibs} = 2$, then the
risk increases two-fold. Finally, a $\lambda_{sibs} $ shows significant as its
95\% confidence interval exceeds 1.  }
\usage{ stats_lambdasibs(data, class) }
\arguments{
   \item{data}{ A data matrix,}
   \item{class}{ A class label.}
}
\author{ Fabrice Colas, Jeanine Houwing-Duistermaat}
