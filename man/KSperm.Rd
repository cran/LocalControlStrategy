\name{KSperm}
\alias{KSperm}
\title{Estimate a Pmax value for the statistical significance of a Kolmogorov-Smirnov D-statistic.}
\description{For a given number of clusters, KSperm() simulates the NULL distribution of LTDs or LRCs
  resulting from Purely Random Clusterings of experimental units within the given data.frame. This NULL
  distribution is discrete due to occurrence of within-cluster TIED values. The observed D-Statistic
  from a previous call to confirm() is then compared with the order statistics computed by KSperm().
  When KSperm() is called immediately after confirm() and the seed value is known, then both the
  Pmax value and the KS-D order statistics of the simulated NULL distribution will be reproducible.   
}
\usage{KSperm(x, reps=100)}
\arguments{
  \item{x}{The output object from a previous call to confirm().}
  \item{reps}{Number of Simulation Replications. Each experimental unit is used at most once within
    each full Replication. No clusters will be empty, but some may be "uninformative".}
}
\details{The numerical value of the Kolmogorov-Smirnov D-statistic from confirm() can be of interest,
  but its "p.value" from ks.test() is badly biased downwards. This bias results because the distribution
  of LTDs or LRCs across clusters is always discrete, due to TIED values within clusters that typically
  also vary in size. Thus, KSperm() generates "reps" additional, independent, NULL values of KS-D and
  computes their order statistics. Finally, KSperm() compares the observed KS-D from confirm() with its
  simulated NULL order statistics to estimate a maximum p.value, denoted by "Pmax" in KSperm() outputs.
}
\value{ An output list object of class KSperm:
  \item{hiclus}{Hierarchical clustering object created by the designated method.}
  \item{dframe}{Name of data.frame containing X, t & Y variables.}
  \item{trtm}{Name of numerical treatment/exposure variable.}
  \item{yvar}{Name of numerical y-Outcome variable.}
  \item{Type}{1 ==> LTDs, otherwise LRCs.}
  \item{reps}{Number of overall Replications, each with the same number, K, of requested clusters.}
  \item{nclus}{Number of clusters requested.}
  \item{units}{Number of experimental units or patients.}
  \item{obsD}{Observed numerical value of KS D-statistic from confirm()}
  \item{Dvec}{Vector of order statistics for simulated NULL KS D-statistics.}
  \item{Pmax}{Estimate of Maximum p.value from discrete LTD/LRC distributions.}
}
\references{
  Obenchain RL. (2010) Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192. 

  Obenchain RL. (2018) \bold{LCstrategy_in_R.pdf} http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{confirm}}, \code{\link{ltdagg}} and \code{\link{lrcagg}}.}
\keyword{nonparametric}