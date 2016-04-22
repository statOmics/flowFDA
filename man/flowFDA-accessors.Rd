\docType{methods}
\name{flowFDA-accessors}
\alias{flowFDA-accessors}
\alias{getBasis}
\alias{getBasis,flowBasis-method}
\alias{getBw}
\alias{getBw,flowBasis-method}
\alias{getBw,flowDa-method}
\alias{getBw,flowPca-method}
\alias{getClust}
\alias{getClust,flowPca-method}
\alias{getClustClass}
\alias{getClustClass,flowPca-method}
\alias{getDa}
\alias{getDa,flowDa-method}
\alias{getDaScore}
\alias{getDaScore,flowDa-method}
\alias{getFset}
\alias{getFset,flowBasis-method}
\alias{getGroups}
\alias{getGroups,flowDa-method}
\alias{getMpc}
\alias{getMpc,flowDa-method}
\alias{getParam}
\alias{getParam,flowBasis-method}
\alias{getParam,flowDa-method}
\alias{getParam,flowPca-method}
\alias{getPbFp}
\alias{getPbFp,flowBasis-method}
\alias{getPbMod}
\alias{getPbMod,flowBasis-method}
\alias{getPca}
\alias{getPca,flowPca-method}
\alias{getPcaScore}
\alias{getPcaScore,flowPca-method}
\alias{nPca}
\alias{nPca,flowDa-method}
\alias{nPca,flowPca-method}
\alias{nPca<-}
\alias{nPca<-,flowPca-method}
\alias{nSet}
\alias{nSet,flowBasis-method}
\alias{nSet,flowDa-method}
\alias{nSet,flowPca-method}
\alias{nbin}
\alias{nbin,flowBasis-method}
\alias{nbin,flowDa-method}
\alias{nbin,flowPca-method}
\alias{setClust}
\alias{setClust<-}
\alias{setClust<-,flowPca-method}
\alias{show}
\title{Extract info from flowBasis, flowPca and flowDa objects}
\arguments{
  \item{x}{a flowBasis, flowPca or flowDa object}
}
\value{
  \code{nbin(x)} returns the number of bins used for each
  flow cytometry channel

  \code{nSet(x)} returns the number of flowSets in the flow
  cytometry experiment

  \code{getParam(x)} returns the flow cytometry channel
  combinations

  \code{getFset(x)} returns the original flowSet with the
  raw flowcytometry data if the flowBasis object was
  generated with flag saveFcs=TRUE or if probability
  binning was used.

  \code{getPbFp(x)} returns the flowFP:flowFP object used
  to setup the basis object in case probability binning was
  used to generate the flowBasis object (using
  probBin=TRUE).

  \code{getPbMod(x)} returns the flowFP:flowFPModel object
  used to setup the flowBasis object when probability
  binning is used to construct the basis (using
  probBin=TRUE).

  \code{getBw(x)} returns the bandwidth of the kernel
  density estimator used to construct the flowBasis object

  \code{getBasis(x)} extracts the basis from the flowBasis
  object

  \code{getPca(x)} extracts the basis from a flowPca or
  flowDa object

  \code{getClust(x)} extracts the clustering object from a
  flowPca or flowDa object

  \code{setClust(x)<-} replace clust of flowPca or flowDa
  object

  \code{getPcaScore(x)} extracts the PCA Scores from a
  flowPca or flowDa object

  \code{getClustClass} extracts the cluster membership from
  the clustering object of a flowPca or flowDa object

  \code{nPca(x)} extracts the number of PCs used for model
  based clustering or constructing the flowDa object

  \code{nPca(x)<-} replaces nPca slot of flowPca or flowFda
  objects

  \code{getMpc(x)} extracts the result of multiple testing
  between groups in the discriminant space

  \code{getDa(x)} extracts the discriminant object from a
  flowDa-class object#'

  \code{getGroups(x)} extracts the grouping

  \code{getDaScore(x)} extracts the discriminant scores
  from a flowDa object
}
\description{
  This file discribes different ways to acces the slots and
  values contained in \code{\link{flowBasis-class}},
  \code{\link{flowPca-class}} and
  \code{\link{flowDa-class}} objects.
}

