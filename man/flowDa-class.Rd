\name{flowDa-class}
\alias{flowDa-class}
\title{The flowDa class}
\description{
  This class extends the flowPca class by constructing a
  Fisher Discriminant Analysis using Fisher's Method
}
\section{Slots}{
  \describe{ \item{\code{daObj}:}{A list with the output
  from the Fisher Discriminant Analysis (FDA)}
  \item{\code{groups}:}{A factor object with the class
  labels for each sample} \item{\code{nPca}:}{An integer
  with the number of principal components that are used in
  the FDA} \item{\code{mpc}:}{A list with multiple
  comparison tests between all classes in the discriminant
  space} \item{\code{pcaObj}:}{A stat:prcomp object with
  output form the principal component analysis using the
  function prcomp} \item{\code{pcom}:}{A matrix with the
  channel combinations for which bivariate densities are
  calculated in the flowBasis object}
  \item{\code{param}:}{A matrix with the names of the
  channel combinations for which bivariate densities are
  calculated in the flowBasis object}
  \item{\code{nbin}:}{An integer indicating the number of
  bins used for the flowBasis object, i.e. an nbin x nbin
  grid is used for constructing the basis}
  \item{\code{bw}:}{The bandwith for the kernel density
  estimator used to construct the flowBasis object}
  \item{\code{probBin}:}{Logical flag to indicate if
  probability binning of flowFP is used for construction of
  basis. Probability binning is provided for compatibility
  with De Roy et al. (2012).} }
}
\references{
  De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N.
  (2012). Flow cytometry for fast microbial community
  fingerprinting. Water Research, 46 (3), 907-919.
}

