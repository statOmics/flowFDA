\name{flowBasis}
\alias{flowBasis}
\title{flowBasis constructor}
\usage{
  flowBasis(fcs, param, nbin = 128, bw = 0.01,
    normalize = function(x) x/max(x), probBin = FALSE,
    saveFcs = FALSE)
}
\arguments{
  \item{fcs}{A flowCore:flowSet-class object for which a
  flowBasis is desired}

  \item{param}{Channel names of flowcytometry experiment
  for which bivariate density basis are derived.}

  \item{nbin}{Number of bins that are taken in each channel
  to approximate the bivariate densities}

  \item{bw}{Bandwidth of for the kernel density estimator
  calculated at each bin}

  \item{normalize}{A user defined function for rescaling or
  normalising the densities, standard the density estimates
  are rescaled between 0 and 1}

  \item{probBin}{for compatibility, Flag to use probability
  binning approach of the flowFP package. Was used in De
  Roy et al. 2012. The probBin approach lacks the
  interpretation feature and requires the basis to be
  recalculated when new samples are added to the
  analysis.#'}

  \item{saveFcs}{logical flag to save original
  flowCore:flowSet object within flowBasis object. The
  default settings are saveFcs=FALSE.(The flowSet object is
  always saved when probability binning is used)}
}
\value{
  An instance of an object of the
  \code{\link{flowBasis-class}}
}
\description{
  Constructor for \code{\link{flowBasis-class}}
}
\examples{
if(require(flowFDAExampleData)){
# load a flowSet to use as an example basis object load(fbasis)
data(fset)
data(group)
param=c("SS Log","FL 1 Log","FL 3 Log")
fbasis<-flowBasis(fset,param)
#flowda<-flowDa(fbasis,group,nPca=6)
#individual densities, average densities can be plotted by using scalar with the sample number or a vector with the sample numbers, respectively
samp=1:6
plot(fbasis,samples=samp)
}
}
\references{
  De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N.
  (2012). Flow cytometry for fast microbial community
  fingerprinting. Water Research, 46 (3), 907-919.
}
\seealso{
  \code{\link{flowBasis-class}}
}

