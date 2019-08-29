\name{mlme}
\alias{mlme}
\title{Create a <<Most-Like-Me>> data.frame for a specified X-Confounder vector: xvec}
\description{For a Given X-confounder Vector (xvec), sort all experimental units (eUnits)
  in an ltdagg() or lrcagg() output object into the strictly non-decreasing order of
  their distances from this X-Vector, which defines the TARGET eUnit: "Me". Plots of
  mlme() objects and displays of mlme.stats() are then used to Visualize and Summarize
  "Mini-" << LOCAL effect-size Distributions >> for different Numbers of "Nearest
  Neighbor" eUnits.
}
\usage{
  mlme(envir, hcl, LCagg, xvec )
}
\arguments{
  \item{envir}{Environment output by a call to the LCsetup() function.}
  \item{hcl}{Name of a LCcluster() output object created using a cluster::diana or
     stats::hclust method.}
  \item{LCagg}{A data.frame object output by ltdagg() or lrcagg() containing LOCAL
     effect-size Estimates for eUnits within Clusters defined in X-covariate space.}
  \item{xvec}{A suitable vector of the Numerical values for the X-Confounder variables,
     used in the current CLUSTERING, that define the eUnit: "Me".}
}
\details{For example, in demo(radon), the eUnits are 2881 US "Counties", and the LCagg
  object is of type lrcagg() because radon exposure is a continuous variable. But, in
  demo(pci15k), the eUnits are 15487 "Patients," and the LCagg object is of type
  ltdagg() because treatment choice (thin) is Binary (0 = "No", 1 = "Yes").
}
\value{ An output list object of class mlme:
  \item{xvec}{The xvec vector input to mlme().}
  \item{Type}{Either "LTD" or "LRC".}
  \item{xvars}{Names of the X-Confounder variables specified in LCsetup().}
  \item{varx}{The vector of Variances of the xvars variabes, used in rescaling distances.} 
  \item{outdf}{The output data.frame of sorted "Nearest Neighbor" candidate eUnits.}
}
\references{
  Obenchain RL. LocalControlStrategy-vignette. (2019) \bold{LCstrategy_in_R.pdf}
  http://localcontrolstatistics.org 
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{plot.mlme}},\code{\link{print.mlme}},\code{\link{mlme.stats}}}
\examples{
  \dontrun{
  # Long running example...
  data(pci15k)
  xvars <- c("stent", "height", "female", "diabetic", "acutemi", "ejfract", "ves1proc")
  hclobj <- LCcluster(pci15k, xvars)
  LC.env <- LCsetup(hclobj, pci15k, thin, surv6mo)
  surv0500 <- ltdagg(500, LC.env)
  xvec11870 <- c( 0, 162, 1, 1, 0, 57, 1) 
  mlmeC5H <- mlme(envir = LC.env, hcl = hclobj, LCagg = surv0500, xvec = xvec11870 ) 	
  plot(mlmeC5H) # using default "NN" and "breaks" settings...
  }
}
\keyword{nonparametric}
