\name{flowBasis-class}
\alias{flowBasis-class}
\title{The flowBasis class}
\description{
  This class represents a flowset by a basis derived from a
  kernel density estimator using an equally spaced grid.
}
\section{Slots}{
  \describe{ \item{\code{basis}:}{A matrix with the basis
  coefficients} \item{\code{pcom}:}{A matrix with the
  channel combinations for which bivariate densities are
  calculated} \item{\code{param}:}{A matrix with the names
  of the channel combinations for which bivariate densities
  are calculated} \item{\code{nbin}:}{An integer indicating
  the number of bins for each variable, i.e. an nbin x nbin
  grid is used for constructing the basis}
  \item{\code{bw}:}{The bandwith that is used for the
  kernel density estimator} \item{\code{probBin}:}{Logical
  flag to indicate if probability binning of flowFP is used
  for construction of basis. Probability binning is
  provided for compatibility with De Roy et al. (2012).}
  \item{\code{fset}:}{save original flowCore:flowSet-class
  object used to construct the basis}
  \item{\code{fmod}:}{flowFPModel used when probBin=TRUE.
  If probBin is FALSE the model is empty and the improved
  approach is adopted.  Probability binning is provided for
  compatibility with De Roy et al. (2012).}
  \item{\code{fp}:}{flowFP fingerprint used when
  probBin=TRUE. If probBin is FALSE the fp object is empty
  and the improved approach is adopted.  Probability
  binning is provided for compatibility with De Roy et al.
  (2012).} }
}
\references{
  De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N.
  (2012). Flow cytometry for fast microbial community
  fingerprinting. Water Research, 46 (3), 907-919.
}

