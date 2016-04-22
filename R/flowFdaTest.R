#' @title Create slot for testing differences between treatments in the discriminant space
#'
#' @aliases flowDaTest
#' @rdname flowDaTest-method
#' @docType methods
#' 
#' @param disc the discriminants for which testing has to be performed
#' @param nPerm the number of permutations that has to be performed for establishing the permutation distributions of differences in the discriminant space
#' @param what="p.value" for establishing the null distribution of the p-values and "what=statistic" for establishing the null distribtion of the statistic.
#'
#' @references De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919. 
#' @rdname flowDaTest
#' @aliases flowDaTest,flowDa-method
setMethod("flowDaTest",signature="flowDa",definition=function(fDa,disc,nPerm=5000,what="p.value")
{
  if(!is.element(what,c("statistic","p.value"))) stop("what-argument should be statistic or p.value")
  ndisc=length(disc)
	###setup tests
	gLevels<-levels(getGroups(fDa))
	nLevels<-nlevels(getGroups(fDa))
	groupMeans<-t(sapply(gLevels,function(x,M,groups) colMeans(M[groups==x,]),groups=getGroups(fDa),M=getDa(fDa)$M))
	mpc<-matrix("",ncol=2,nrow=nLevels*(nLevels-1)/2)
	teller=0
	for (i in 1:(nLevels-1)) 
		for (j in (i+1):(nLevels)) 
			{ 
				teller=teller+1
				mpc[teller,]<-gLevels[c(i,j)]	
			}
	pOrig<-.mpcAll(getDa(fDa)$M,getGroups(fDa),mpc,disc,what)
	###setup permutated linear model
  nProgBlocks<-nPerm%/%500
  cat("progress",rep("*",nProgBlocks),"\n        ")
	pPerm=sapply(1:nPerm,function(x,fDa,disc,mpc,what) 
	{
  if(x%%500==0) cat(" *")
	id<-sample(nrow(getDa(fDa)$M))
        ds=ncol(getDa(fDa)$M)
	Lm<-lm(getPca(fDa)$x[id,1:nPca(fDa)]~-1+getGroups(fDa))
	W<-t(getPca(fDa)$x[id,1:nPca(fDa)]-Lm$fitted)%*%(getPca(fDa)$x[id,1:nPca(fDa)]-Lm$fitted)/(length(id)-nlevels(getGroups(fDa)))
	svdQInv<-svd(diag(1/getPca(fDa)$sdev[1:nPca(fDa)]^2)%*%W) 
	M=(getPca(fDa)$x[id,1:nPca(fDa)]%*%svdQInv$v[,nPca(fDa):1])[,1:ds]
	return((.mpcAll(M,getGroups(fDa),mpc,disc,what)))
	},fDa=fDa,disc=disc,mpc=mpc,what=what)
	
	#calculate permulated pvalue
	if (what=="p.value") pValuePerm<-matrix(rowMeans(pPerm<c(pOrig)),ncol=length(disc))
  else pValuePerm<-matrix(rowMeans(abs(pPerm)>abs(c(pOrig))),ncol=length(disc))
	rownames(pValuePerm)<-paste(mpc[,2],mpc[,1],sep="-")
	colnames(pValuePerm)<-paste("D",disc,sep="")
	#idem for contrast matrix L
	contrasts=contrMat(1:nLevels,type="Tukey")%*%groupMeans[,disc]
	colnames(contrasts)<-paste("D",disc,sep="")
	fDa@mpc<-list(contrasts=contrasts,pValuePerm=pValuePerm)	
  cat("\n")
return(fDa)
})


