#########accessors and show methods
#'Extract info from flowBasis, flowPca and flowDa objects
#'
#' This file discribes different ways to acces the slots and values contained in \code{\link{flowBasis-class}}, \code{\link{flowPca-class}} and  \code{\link{flowDa-class}} objects.
#'
#' @param x a flowBasis, flowPca or flowDa object 
#'
#' @name flowFDA-accessors
#' @rdname flowFDA-accessors
#' @aliases nbin getParam getBw getBasis getPca show getClust nPca getMpc getGroups getDa getClustClass setClust<- nPca<- getPbMod getFset getPbFp
#' @docType methods


#' @return \code{nbin(x)} returns the number of bins used for each flow cytometry channel
#' @rdname flowFDA-accessors
#' @aliases nbin,flowBasis-method
setMethod("nbin",signature="flowBasis",definition=function(x) return(x@nbin))
#' @rdname flowFDA-accessors
#' @aliases nbin,flowPca-method
setMethod("nbin",signature="flowPca",definition=function(x) return(x@nbin))
#' @rdname flowFDA-accessors
#' @aliases nbin,flowDa-method
setMethod("nbin",signature="flowDa",definition=function(x) return(x@nbin))


#' @return \code{nSet(x)} returns the number of flowSets in the flow cytometry experiment
#' @rdname flowFDA-accessors
#' @aliases nSet,flowBasis-method
setMethod("nSet",signature="flowBasis",definition=function(x) return(nrow(x@basis)))
#' @rdname flowFDA-accessors
#' @aliases nSet,flowPca-method
setMethod("nSet",signature="flowPca",definition=function(x) return(nrow(x@pcaObj$x)))
#' @rdname flowFDA-accessors
#' @aliases nSet,flowDa-method
setMethod("nSet",signature="flowDa",definition=function(x) return(nrow(x@pcaObj$x)))

#' @return \code{getParam(x)} returns the flow cytometry channel combinations
#' @rdname flowFDA-accessors
#' @aliases getParam,flowBasis-method
setMethod("getParam",signature="flowBasis",definition=function(x) return(x@param))
#' @rdname flowFDA-accessors
#' @aliases getParam,flowPca-method
setMethod("getParam",signature="flowPca",definition=function(x) return(x@param))
#' @rdname flowFDA-accessors
#' @aliases getParam,flowDa-method
setMethod("getParam",signature="flowDa",definition=function(x) return(x@param))

#'@return \code{getFset(x)} returns the original flowSet with the raw flowcytometry data if the flowBasis object was generated with flag saveFcs=TRUE or if probability binning was used.  
#'@rdname flowFDA-accessors
#'@aliases getFset,flowBasis-method
setMethod("getFset",signature="flowBasis",definition=function(x) return(x@fset))

#'@return \code{getPbFp(x)} returns the flowFP:flowFP object used to setup the basis object in case probability binning was used to generate the flowBasis object (using probBin=TRUE).  
#'@rdname flowFDA-accessors
#'@aliases getPbFp,flowBasis-method
setMethod("getPbFp",signature="flowBasis",definition=function(x) return(x@fp))

#'@return \code{getPbMod(x)} returns the flowFP:flowFPModel object used to setup the flowBasis object when probability binning is used to construct the basis (using probBin=TRUE).  
#'@rdname flowFDA-accessors
#'@aliases getPbMod,flowBasis-method
setMethod("getPbMod",signature="flowBasis",definition=function(x) return(x@fmod))


#' @return \code{getBw(x)} returns the bandwidth of the kernel density estimator used to construct the flowBasis object
#' @rdname flowFDA-accessors
#' @aliases getBw,flowBasis-method
setMethod("getBw",signature="flowBasis",definition=function(x) return(x@bw))
#' @rdname flowFDA-accessors
#' @aliases getBw,flowPca-method
setMethod("getBw",signature="flowPca",definition=function(x) return(x@bw))
#' @rdname flowFDA-accessors
#' @aliases getBw,flowDa-method
setMethod("getBw",signature="flowDa",definition=function(x) return(x@bw))

#' @return \code{getBasis(x)} extracts the basis from the flowBasis object
#' @rdname flowFDA-accessors
#' @aliases getBasis,flowBasis-method
setMethod("getBasis",signature="flowBasis",definition=function(x) return(x@basis))

#' @return \code{getPca(x)} extracts the basis from a flowPca or flowDa object
#' @rdname flowFDA-accessors
#' @aliases getPca,flowPca-method
setMethod("getPca",signature="flowPca",definition=function(x) return(x@pcaObj))

#' @return \code{getClust(x)} extracts the clustering object from a flowPca or flowDa object
#' @rdname flowFDA-accessors
#' @aliases getClust,flowPca-method
setMethod("getClust",signature="flowPca",definition=function(x) return(x@clust))
#' @return \code{setClust(x)<-} replace clust of flowPca or flowDa object
#' @name setClust
#' @aliases setClust<-,flowPca-method
#' @rdname flowFDA-accessors
setReplaceMethod("setClust",signature(x="flowPca"),function(x,value){
x@clust<-value
validObject(x)
x
})


#' @return \code{getPcaScore(x)} extracts the PCA Scores from a flowPca or flowDa object
#' @rdname flowFDA-accessors
#' @aliases getPcaScore,flowPca-method
setMethod("getPcaScore",signature="flowPca",definition=function(x,nPca=NULL) if(!is.null(nPca)) {if( (nPca>0)&(ncol(x@pcaObj$x)>=nPca)) return(x@pcaObj$x[,1:nPca]) else return(x@pcaObj$x)} else return(x@pcaObj$x))

#' @return \code{getClustClass} extracts the cluster membership from the clustering object of a flowPca or flowDa object
#' @rdname flowFDA-accessors
#' @aliases getClustClass,flowPca-method
setMethod("getClustClass",signature="flowPca",definition=function(x) return(x@clust$classification))


#' @return \code{nPca(x)} extracts the number of PCs used for model based clustering or constructing the flowDa object
#' @rdname flowFDA-accessors
#' @aliases nPca,flowPca-method
setMethod("nPca",signature="flowPca",definition=function(x) return(x@nPca))

#' @return \code{nPca(x)<-} replaces nPca slot of flowPca or flowFda objects
#' @name nPca
#' @aliases nPca<-,flowPca-method
#' @rdname flowFDA-accessors

setReplaceMethod("nPca",signature(x="flowPca"),definition=function(x,value) {
if (value<1) x@nPca=as.numeric(which((cumsum(x@pcaObj$sdev^2)/sum(x@pcaObj$sdev^2))>value)[1]) else (x@nPca=value)
if (x@nPca<2) x@nPca=2
validObject(x)
x
}) 
#' @rdname flowFDA-accessors
#' @aliases nPca,flowDa-method
setMethod("nPca",signature="flowDa",definition=function(x) return(x@nPca))

#' @return \code{getMpc(x)} extracts the result of multiple testing between groups in the discriminant space
#' @rdname flowFDA-accessors
#' @aliases getMpc,flowDa-method 
setMethod("getMpc",signature="flowDa",definition=function(x) return(x@mpc))

#' @return \code{getDa(x)} extracts the discriminant object from a flowDa-class object#' @rdname flowFDA-accessors
#' @aliases getDa,flowDa-method 
setMethod("getDa",signature="flowDa",definition=function(x) return(x@daObj))

#' @return \code{getGroups(x)} extracts the grouping 
#' @rdname flowFDA-accessors
#' @aliases getGroups,flowDa-method 
setMethod("getGroups",signature="flowDa",definition=function(x) return(x@groups))


#' @return \code{getDaScore(x)} extracts the discriminant scores from a flowDa object
#' @rdname flowFDA-accessors
#' @aliases getDaScore,flowDa-method
setMethod("getDaScore",signature="flowDa",definition=function(x) return(getDa(x)$M))


