#####Setup Basis using kernel smoother
.basisSetup<-function(set,nbin,bw)
{
  return(grDevices:::.smoothScatterCalcDensity(rbind(set,c(1,1),c(0,0)),nbin,bw)$fhat)
}

#' flowBasis constructor
#' 
#' Constructor for \code{\link{flowBasis-class}}
#' 
#' @param fcs A flowCore:flowSet-class object for which a  flowBasis is desired
#' @param param Channel names of flowcytometry experiment for which bivariate density basis are derived.
#' @param nbin Number of bins that are taken in each channel to approximate the bivariate densities
#' @param bw Bandwidth of for the kernel density estimator calculated at each bin
#' @param normalize A user defined function for rescaling or normalising the densities, standard the density estimates are rescaled between 0 and 1
#' @param probBin for compatibility, Flag to use probability binning approach of the flowFP package. Was used in De Roy et al. 2012. The probBin approach lacks the interpretation feature and requires the basis to be recalculated when new samples are added to the analysis.#' @param saveFcs logical flag to save original flowCore:flowSet object within flowBasis object. The default settings are saveFcs=FALSE.(The flowSet object is always saved when probability binning is used)    
#' @return An instance of an object of the \code{\link{flowBasis-class}}
#' @examples 
#' if(require(flowFDAExampleData)){ 
#' # load a flowSet to use as an example basis object load(fbasis)
#' data(fset)
#' data(group)
#' param=c("SS Log","FL 1 Log","FL 3 Log")
#' fbasis<-flowBasis(fset,param)
#' #flowda<-flowDa(fbasis,group,nPca=6)
#' #individual densities, average densities can be plotted by using scalar with the sample number or a vector with the sample numbers, respectively
#' samp=1:6
#' plot(fbasis,samples=samp)
#' }
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#' @seealso \code{\link{flowBasis-class}}
#' @importFrom flowCore 
#' @export
flowBasis<-function(fcs,param,nbin=128,bw=0.01,normalize=function(x) x/max(x),probBin=FALSE,saveFcs=FALSE)
{
fBasis<-new("flowBasis")
p<-length(param)
fBasis@pcom<-matrix(0,ncol=2,nrow=(p-1)*p/2)
teller<-0
for (i in 1:(p-1))
for (j in (i+1):p)
{
teller<-teller+1
fBasis@pcom[teller,]<-c(i,j)
}
if (probBin)
{
fBasis@param=matrix(param[fBasis@pcom],ncol=2)
fBasis@nbin=as.integer(nbin)
fBasis@fmod=flowFPModel(fcs,parameters=param,nRecursions=ceiling(log2(nbin)))
fBasis@probBin=probBin
fBasis@fp=flowFP(fcs,fBasis@fmod)
fBasis@basis=counts(fBasis@fp)
fBasis@basis=counts(fBasis@fp)/rowSums(fBasis@basis)
fBasis@fset=fcs
return(fBasis)
} else
{
fBasis@probBin=probBin
fBasis@basis=fsApply(fcs,function(x,param,pcom,nbin) {c(sapply(1:nrow(pcom), function(k,set,param,pcom,nbin,bw) c(.basisSetup(set[,param[pcom[k,]]],nbin,bw)),set=exprs(x),param=param,pcom=fBasis@pcom,nbin=nbin,bw=bw))},param=param,pcom=fBasis@pcom,nbin=nbin)
fBasis@param=matrix(param[fBasis@pcom],ncol=2)
fBasis@nbin=as.integer(nbin)
fBasis@bw=bw
fBasis@basis<-normalize(fBasis@basis)
if(saveFcs) fBasis@fset=fcs
return(fBasis)}
}

#' flowPca constructor
#' 
#' Constructor for \code{\link{flowPca-class}}
#' 
#' A flowPca object is a reduced representation of a \code{\link{flowBasis-class}} and a flowCore:flowSet-class with respect to the bivariate distributions of the flow channels that are considered. 
#' Consists of a principal component decomposition of the coefficients of the estimated bivariate kernel densities.
#' The PCs can be interpreted with respect to the features in the original bivariate distributions.
#' 
#' @param fbasis A \code{\link{flowBasis-class}} object for which a discriminant basis is desired.
#' 
#' @return An instance of an object of the \code{\link{flowPca-class}}
#' @seealso \code{\link{flowBasis-class}}, \code{\link{flowPca-class}} and \code{\link{flowDa-class}}
#'
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#'
#' @examples
#' # load a flowSet to use as an example basis object
#' if(require(flowFDAExampleData)){
#' data(fbasis)
#' data(group)
#' data(param)
#' 
#' flowpca<-flowPca(fbasis)
#' plot(flowpca)
#'
#' ####Plot of the loadings the first PC of the daobject
#' pc=1
#' plot(flowpca,disc=pc,main=paste("PC",pc),plotType="pcaInt")
#' 
#' #interpretation of score on PC1 for first flowcytometry sample
#' #pc=1
#' #L=rep(0,length(group))
#' #L[1]<-1
#' #plot(flowpca,flowbasis,disc=pc,L=L,plotType="pcaCont")
#' }
#' @importFrom flowCore 
#' @export

flowPca<-function(fbasis)
{
fPca<-new("flowPca")
fPca@pcaObj<-prcomp(fbasis@basis)
fPca@param=fbasis@param
fPca@pcom=fbasis@pcom
fPca@nbin=fbasis@nbin
fPca@bw=fPca@bw
fPca@probBin=fbasis@probBin
return(fPca)
}

#' flowDa constructor
#' 
#' Constructor for \code{\link{flowDa-class}}
#'
#' A flowDa object is a reduced representation of a \code{\link{flowBasis-class}} and a flowCore:flowSet-class with respect to the bivariate distributions of the flow channels that are considered. 
#' The discriminants are derived by adopting Fisher's Method on the first \code{nPca} principal components that were derived from the coefficients of the bivariate kernel density estimation. 
#' The discriminants can be interpreted with respect to the features in the original bivariate distributions. 
#' 
#' @param fbasis A \code{\link{flowBasis-class}} object for which a discriminant basis is desired.
#' @param groups A vector of the numeric or factor type with the class labels
#' @param nPca The number of PCs or the fraction of the total variation to be used in the PCA dimension reduction prior to the discrimination analysis
#' 
#' @return An instance of an object of the \code{\link{flowDa-class}}
#' 
#' @seealso \code{\link{flowBasis-class}}, \code{\link{flowPca-class}}, \code{\link{flowDa-class}}, \code{\link{flowBasis}} and \code{\link{flowPca}}
#' 
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#'
#' @examples 
#' if(require(flowFDAExampleData)){
#' # load a flowSet to use as an example basis object load(fbasis)
#' data(fbasis)
#' data(group)
#' data(param)
#' flowda<-flowDa(fbasis,group,nPca=6)
#' plot(flowda)
#' 
#' ####Plot of the loadings the first discriminant of the daobject
#' #discriminant=1
#' #plot(flowda,disc=discriminant,plotType="discInt",main=paste("Discriminant",discriminant),cex.axis=1.5,cex.lab=2,cex.main=2)
#' 
#' ####interpretation plot using a vector L to assess the contribution of centroid of the group 1 on the first discrimant
#' ####Other definitions of L support the interpretation of the score for an individual sample or for contrasts between groups etc.
#' #discriminant=1
#' #L=group==levels(group)[1]
#' #plot(flowda,fbasis,disc=discriminant,plotType="discCont",L=L,main=paste("Discriminant",discriminant),cex.axis=1.5,cex.lab=2,cex.main=2)
#' }
#' @importFrom flowCore 
#' @export

flowDa<-function(fbasis,groups,nPca=2)
{
X=model.matrix(~-1+groups) 
fDa<-new("flowDa")
####PCA
fDa@pcaObj<-prcomp(fbasis@basis)
if (nPca<1) nPca=which((cumsum(fDa@pcaObj$sdev^2)/sum(fDa@pcaObj$sdev^2))>nPca)[1]
if (nPca<2) nPca=2
####Fisher's Method
#svd of \hat{\mb{\Sigma}}_within^{-1} \mb{B} =  \hat{\mb{\Sigma}}_within^{-1} \mb{B}
#same solution if you replace B with varcovar matrix Sigma^2_{total}. 
#PCA that is the diag matrix of the eigenvalues
#Q<-ginv(W) %*% diag(fDa@pcaObj$stdev^2)
Lm<-lm(fDa@pcaObj$x[,1:nPca]~-1+X)
W<-t(fDa@pcaObj$x[,1:nPca]-Lm$fitted)%*%(fDa@pcaObj$x[,1:nPca]-Lm$fitted)/(nrow(X)-ncol(X))
#if (ridgePen) svdQ<-svd(solve(W+diag(ncol(W))*ridgePen) %*% diag(fDa@pcaObj$sdev[1:nPca]^2))
#else 
#svdQ<-svd(solve(W) %*% diag(fDa@pcaObj$sdev[1:nPca]^2))
{svdQInv<-svd(diag(1/fDa@pcaObj$sdev[1:nPca]^2)%*%W) #(u[,nPca:1] is v and v[,nPca:1] is u,  1/d[nPca:1] is d)
#svdQ<-svdQInv
#svdQ$d<-1/svdQInv$d[nPca:1]
#svdQ$u<-svdQInv$v[,nPca:1]
#svdQ$v<-svdQInv$u[,nPca:1]
}
ds=min(nrow(X),ncol(X)-1)
delta<-svdQInv$v[,nPca:1]
lambda=1/svdQInv$d[nPca:1]
fDa@daObj<-list(delta=delta[,1:ds],lambda=lambda[1:ds],M=fDa@pcaObj$x[,1:nPca]%*%delta[,1:ds])

####Other output
fDa@groups=groups
fDa@nPca=as.integer(nPca)
fDa@param=fbasis@param
fDa@pcom=fbasis@pcom
fDa@nbin=fbasis@nbin
fDa@bw=fbasis@bw
fDa@probBin=fbasis@probBin
return(fDa)
}
