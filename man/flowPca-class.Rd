\name{flowPca-class}
\alias{flowPca-class}
\title{The flowPca class}
\description{
  This class represents a principal component basis derived
  from a flowBasis-class object
}
\section{Slots}{
  \describe{ \item{\code{pcaObj}:}{A stats::prcomp object
  with output form the principal component analysis using
  the function prcomp} \item{\code{pcom}:}{A matrix with
  the channel combinations for which bivariate densities
  are calculated in the flowBasis object}
  \item{\code{param}:}{A matrix with the names of the
  channel combinations for which bivariate densities are
  calculated in the flowBasis object}
  \item{\code{nbin}:}{An integer indicating the number of
  bins used for the flowBasis object, i.e. an nbin x nbin
  grid is used for constructing the basis}
  \item{\code{bw}:}{The bandwith for the kernel density
  estimator used to construct the flowBasis object}
  \item{\code{clust}:}{An mclust:Mclust object with output
  of model-based clustering} \item{\code{nPca}:}{An integer
  with the number of principal components that are used in
  model based clustering, empty if no clustering has been
  done yet} \item{\code{probBin}:}{Logical flag to indicate
  if probability binning of flowFP is used for construction
  of basis. Probability binning is provided for
  compatibility with De Roy et al. (2012).} }
}
\references{
  De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N.
  (2012). Flow cytometry for fast microbial community
  fingerprinting. Water Research, 46 (3), 907-919.
}

