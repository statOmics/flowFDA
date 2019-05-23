###
#use S3 classes prcomp and Mclust as S4 slots 
#prcomp constructor
hlp<-rnorm(100)
dim(hlp)<-c(10,10)

prcompProt=prcomp(hlp)
for (i in (1:length(prcompProt)))
prcompProt[[i]]=new(Class=class(prcompProt[[i]]))

mclustProt=Mclust(hlp)
for (i in (1:length(mclustProt)))
{
className=class(mclustProt[[i]])
mclustProt[[i]]=try(new(Class=class(mclustProt[[i]])),silent=TRUE)
if (class(mclustProt[[i]])=="try-error") 
{	
#mclustProt[[i]]=NULL	
mclustProt[[i]]= "not initialized"
class(mclustProt[[i]])=className
}
}

setOldClass("prcomp",prototype=prcompProt)
setOldClass("Mclust",prototype=mclustProt)
rm(hlp,prcompProt,mclustProt)


# S4 documentation using Roxygen.
#' The flowBasis class
#' 
#' This class represents a flowset by a basis derived from a kernel density estimator using an equally spaced grid.
#'
#' @section Slots: 
#' \describe{
#'  \item{\code{basis}:}{A matrix with the basis coefficients}
#'  \item{\code{pcom}:}{A matrix with the channel combinations for which bivariate densities are calculated}
#'  \item{\code{param}:}{A matrix with the names of the channel combinations for which bivariate densities are calculated}
#' \item{\code{nbin}:}{An integer indicating the number of bins for each variable, i.e. an nbin x nbin grid is used for constructing the basis}
#' \item{\code{bw}:}{The bandwith that is used for the kernel density estimator}
#' \item{\code{probBin}:}{Logical flag to indicate if probability binning of flowFP is used for construction of basis. Probability binning is provided for compatibility with De Roy et al. (2012).}
#' \item{\code{fset}:}{save original flowCore:flowSet-class object used to construct the basis}
#' \item{\code{fmod}:}{flowFPModel used when probBin=TRUE. If probBin is FALSE the model is empty and the improved approach is adopted.  Probability binning is provided for compatibility with De Roy et al. (2012).}
#' \item{\code{fp}:}{flowFP fingerprint used when probBin=TRUE. If probBin is FALSE the fp object is empty and the improved approach is adopted.  Probability binning is provided for compatibility with De Roy et al. (2012).}
#' }
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#' @name flowBasis-class
#' @rdname flowBasis-class
#' @exportClass flowBasis
#' 

setClass("flowBasis",representation(basis="matrix",param="matrix",pcom="matrix",nbin="integer",bw="numeric",probBin="logical",fset="flowSet",fmod="flowFPModel",fp="flowFP"),prototype = list(basis=matrix(nrow=0,ncol=0),param=matrix(ncol=2,nrow=0),pcom=matrix(ncol=2,nrow=0),nbin=integer(),bw=numeric(),probBin=logical(),fset=new(Class="flowSet"),fmod=new(Class="flowFPModel"),fp=new(Class="flowFP")))

#' The flowPca class
#' 
#' This class represents a principal component basis derived from a flowBasis-class object
#'
#' @section Slots: 
#'  \describe{
#'  \item{\code{pcaObj}:}{A stats::prcomp object with output form the principal component analysis using the function prcomp}
#'  \item{\code{pcom}:}{A matrix with the channel combinations for which bivariate densities are calculated in the flowBasis object}
#'  \item{\code{param}:}{A matrix with the names of the channel combinations for which bivariate densities are calculated in the flowBasis object}
#' \item{\code{nbin}:}{An integer indicating the number of bins used for the flowBasis object, i.e. an nbin x nbin grid is used for constructing the basis}
#' \item{\code{bw}:}{The bandwith for the kernel density estimator used to construct the flowBasis object}
#' \item{\code{clust}:}{An mclust:Mclust object with output of model-based clustering}
#' \item{\code{nPca}:}{An integer with the number of principal components that are used in model based clustering, empty if no clustering has been done yet}
#' \item{\code{probBin}:}{Logical flag to indicate if probability binning of flowFP is used for construction of basis. Probability binning is provided for compatibility with De Roy et al. (2012).}
#' }
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#' @name flowPca-class
#' @rdname flowPca-class
#' @exportClass flowPca
setClass("flowPca",representation(pcaObj="prcomp",param="matrix",pcom="matrix",nbin="integer",bw="numeric",clust="Mclust",nPca="numeric",probBin="logical"),prototype=list(pcaObj=new(Class="prcomp"),param=matrix(ncol=2,nrow=0),pcom=matrix(ncol=2,nrow=0),nbin=integer(),bw=numeric(),clust=new(Class="Mclust"),nPca=numeric(),probBin=logical()))

#' The flowDa class
#' 
#' This class extends the flowPca class by constructing a Fisher Discriminant Analysis using Fisher's Method
#' 
#' @section Slots: 
#'  \describe{
#'  \item{\code{daObj}:}{A list with the output from the Fisher Discriminant Analysis (FDA)}
#'  \item{\code{groups}:}{A factor object with the class labels for each sample}
#'  \item{\code{nPca}:}{An integer with the number of principal components that are used in the FDA}
#'  \item{\code{mpc}:}{A list with multiple comparison tests between all classes in the discriminant space}
#'  \item{\code{pcaObj}:}{A stat:prcomp object with output form the principal component analysis using the function prcomp}
#'  \item{\code{pcom}:}{A matrix with the channel combinations for which bivariate densities are calculated in the flowBasis object}
#'  \item{\code{param}:}{A matrix with the names of the channel combinations for which bivariate densities are calculated in the flowBasis object}
#' \item{\code{nbin}:}{An integer indicating the number of bins used for the flowBasis object, i.e. an nbin x nbin grid is used for constructing the basis}
#' \item{\code{bw}:}{The bandwith for the kernel density estimator used to construct the flowBasis object}
#' \item{\code{probBin}:}{Logical flag to indicate if probability binning of flowFP is used for construction of basis. Probability binning is provided for compatibility with De Roy et al. (2012).}
#' }
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#' @name flowDa-class
#' @rdname flowDa-class
#' @exportClass flowDa
setClass("flowDa",representation(daObj="list",groups="factor",mpc="list"),contains="flowPca",prototype=list(daObj=list(),groups=factor(),mpc=list()))

