\name{confirm}
\alias{confirm}
\title{Confirm that Clustering in Covariate X-space yields an "adjusted" LTD/LRC effect-size Distribution}
\description{For a given Number of Clusters, K, confirm() compares the observed distribution of LTDs or
  LRCs from relatively well-matched experimental units with the corresponding distribution from Purely
  Random Clusterings of experimental units. The larger are differences between the (blue) observed
  empirical CDF of effect-sizes and the (red) Purely Random CDF, the more potentially IMPORTANT are the
  "adjustments" resulting from focussing upon clustering (matching) of experimental units in X-space.   
}
\usage{confirm(x, reps=100, seed=12345)}
\arguments{
  \item{x}{An output object from ltdagg() or lrcagg() for a specified number of clusters, K.}
  \item{reps}{Number of simulation Replications, each with the same number, K, and sizes, N1,
    N2, ..., NK of Purely Random clusters.}
  \item{seed}{This (arbitrary) integer argument will be passed to the R set.seed() function.
    Knowing the value of this seed makes the output from confirm() reproducible.}
}
\details{Making multiple calls to confirm() for ltdagg() or lrcagg() objects resulting from
  different choices of K = Numbers of Clusters help the analyst decide which observed LTD or LRC
  effect-size distributions are (or are not) meaningfully different from Purely Random. When the
  X-covariates used in LCcluster() are truly "ignorable," then [i] the resulting clusters will
  be Purely Random, and [ii] both the number (K) and the sizes (N1, N2, ...,NK) of clusters
  formed will be meaningless and arbitrary. Thus confirm() simulates the empirical CDF for LTDs or
  LRCs resulting from purely random permutations of the Cluster ID numbers (1, 2, ...,K) assigned
  by ltdagg() or lrcagg(). Each permutation yields K artificial "clusters" of sizes N1, N2, ..., NK.
  Simulation results are accumulated for the total number of random permutations specified in the
  "reps=" argument of confirm(). Calls to print.confirm() and plot.confirm() provide information on
  comparisons of empirical CDFs for the observed and Purely Random LTD/LRD distributions, including
  output from ks.test(). While the value of the Kolmogorov-Smirnov D-statistic can be of interest,
  the "p.value" from ks.test() is badly biased downwards because LTD/LRC distributions contain MANY
  TIES. Analysts interested in simulating a crude "pmax.value" can invoke KSperm(confirm()), but
  each such follow-up computation is slightly more intensive than the original confirm() computation.
}
\value{ An output list object of class confirm:
  \item{hiclus}{Hierarchical clustering object created by the designated method.}
  \item{dframe}{Name of data.frame containing X, trex & Y variables.}
  \item{trtm}{Name of numerical trex variable.}
  \item{yvar}{Name of numerical Y-outcome variable.}
  \item{reps}{Number of overall Replications, each with the same numbers of requested clusters.}
  \item{seed}{Integer argument passed to set.seed(). Knowing which seed value was used in the
      call to confirm() makes not only the NULL distribution of observed LTDs or LRCs reproducible
	  but also makes the NULL distribution of D-statistics (adjusted for ties) from a subsequent
	  call to KSperm() reproducible.}
  \item{nclus}{Number of clusters requested.}
  \item{units}{Number of experimental units or patients.} 
  \item{Type}{1 ==> LTDs, otherwise LRCs.}
  \item{LCmean}{Weighted Local Mean across Clusters.} 
  \item{LCstde}{Weighted Std. Error across Clusters.} 
  \item{RPmean}{Weighted Random Permutation Mean across Clusters.} 
  \item{RPstde}{Weighted Random Permutation Std. Error across Clusters.}  
  \item{KSobsD}{Output from print(ks.test()).}
  \item{LCdist}{data.frame of 5 key variables for all experimental units.}
  \item{dfconf}{data.frame of lstat = LTD or LRC values of max(length) = reps*units.}
}
\references{
  Obenchain RL. (2010) The Local Control Approach using JMP. Chapter 7 of
  \bold{Analysis of Observational Health Care Data using SAS}, \emph{Cary, NC:SAS Press},
  pages 151-192. 

  Obenchain RL. (2018) \bold{LCstrategy_in_R.pdf} http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{ltdagg}} and \code{\link{lrcagg}}.}
\keyword{nonparametric}