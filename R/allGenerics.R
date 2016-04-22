#' @export
setGeneric("nbin", function(x) standardGeneric("nbin"))

#' @export
setGeneric("nSet",function(x) standardGeneric("nSet"))

#' @export
setGeneric("getParam", function(x) standardGeneric("getParam"))

#' @export
setGeneric("getBw", function(x) standardGeneric("getBw"))

#' @export
setGeneric("getBasis",function(x) standardGeneric("getBasis"))

#' @export
setGeneric("getPca", function(x) standardGeneric("getPca"))
#' @export
setGeneric("getPcaScore", function(x,...) standardGeneric("getPcaScore"))
#' @export
setGeneric("getDaScore", function(x,...) standardGeneric("getDaScore"))



#' @export
setGeneric("nPca", function(x) standardGeneric("nPca"))
#' @export
setGeneric("nPca<-", function(x,value) standardGeneric("nPca<-"))

#' @export
setGeneric("getMpc", function(x) standardGeneric("getMpc"))

#' @export
setGeneric("getDa", function(x) standardGeneric("getDa"))

#' @export
setGeneric("getFset", function(x) standardGeneric("getFset"))

#' @export
setGeneric("getPbMod", function(x) standardGeneric("getPbMod"))

#' @export
setGeneric("getPbFp", function(x) standardGeneric("getPbFp"))


#' @export
setGeneric("getGroups", function(x) standardGeneric("getGroups"))

setGeneric("flowDaCv", function(fset,param,cv,groups,nbins,pcs,bw,file) standardGeneric("flowDaCv"))

#' @export
setGeneric("setClust<-", function(x,value) standardGeneric("setClust<-"))

#' @export
setGeneric("getClust",function(x) standardGeneric("getClust"))

#' @export
setGeneric("flowDaTest", function(fDa,...) standardGeneric("flowDaTest"))

#' @export
setGeneric("getClustClass",function(x,...) standardGeneric("getClustClass"))
