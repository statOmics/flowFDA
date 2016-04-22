
#' Plot methods
#' Generate plots for  \code{\link{flowBasis-class}}, \code{\link{flowPca-class}} or  \code{\link{flowDa-class}} objects object
#'
#' @aliases plot
#' @rdname plot-methods
#' @docType methods
#' 
#' @param x \code{\link{flowBasis-class}}, \code{\link{flowPca-class}} or  \code{\link{flowDa-class}} objects object.
#' @param y obsolete used for compatibility with S3 plot function.
#' @param samples the flows that will be used for constructing the two dimensional projections of the multivariate distribution. samples=1:5 will plot the average projection for the first 5 samples. To be used for flowBasis method.
#' @param pow parameter for the coloring cfr. the function smoothScatter.
#' @param ask if TRUE the user will be asked to press a button before the next plot is displayed.
#' @param L a matrix with a contrasts between flows for constructing contrast plots. The different contrasts are organized in the columns. 
#' @param contour can be used to add contours on the plots. The default is FALSE.
#' @param main=main,... additional arguments can be passed to customize the plots. e.g. cex, cex.axis, etc. 
#' @return Noting. Side-effect: plot graphs.
#' @examples
#' ###########################################
#' ###flowBasis plots
#' ###########################################
#'
#'
#' if(require(flowFDAExampleData)){
#' data(fbasis)
#' #plot of first 6 samples
#' plot(fbasis,samples=1:6)
#'
#' #Contrast between control samples and nutrient 3h
#' L=rep(0,30)
#' L[1:6]=-1/6
#' L[25:30]=1/6
#' plot(fbasis,L=L)
#' }
#' @rdname plot-methods
#' @aliases plot,flowBasis-method

setMethod("plot","flowBasis",function(x,samples=1,pow=.75,ask=TRUE,L=NULL,colorLim=NULL,colRamp=NULL,contour=FALSE,main=NULL,set=NULL,...)
{
  if (is.null(L)) {
    if (x@probBin){
    if (is.logical(samples)) samples=which(samples)
    if (length(samples)>1) hlpOrig<-colMeans(getBasis(x)[samples,]) else hlpOrig<-(getBasis(x)[samples,])
    if (is.null(set)) set=samples 
    hlp<-hlpOrig-mean(hlpOrig)
    if (is.null(colorLim)) colorLim=max(abs(hlp))
    if  (is.null(colRamp)) colRamp=.blueWhiteRed(41)
     hlpCol=.binColsCont(hlp,colorLim,colRamp,plotOrderAbs=FALSE)
     .myFpDotPlot(x@fp,x@fmod,x@fset,hlpCol$binCol,hlpCol$plotOrder,set=set,param=param,ask=ask,main=main)
     .myFpContPlot(c(hlpOrig),hlpCol$binCol,type="l",main=main,ylab="Density in bin",...)
}
    else{
    if (length(samples)>1) hlp<-colMeans(getBasis(x)[samples,]) else hlp<-(getBasis(x)[samples,]) 
    par(ask=ask)
    if (is.null(colRamp)) colRamp=.jet.colors(512)
    for (i in 1:nrow(getParam(x)))
    {
      hlp2<-matrix(hlp[1:(nbin(x)*nbin(x))+(i-1)*(nbin(x)*nbin(x))],ncol=nbin(x))
      image(hlp2^pow,col=colRamp,xlab=getParam(x)[i,1],ylab=getParam(x)[i,2],main=main,...)
    }
    }
  } else
  {
    if(class(L) !="matrix") L<-matrix(L)
    if(is.null(x)) stop("x object is not defined")
    if (nrow(L)!=nrow(getBasis(x))) stop("contrast L has the wrong dimension")
    if (is.null(colorLim))
	{ for (k in 1:ncol(L))
          colorLim=max(abs(L[,k]%*%getBasis(x)),colorLim)
        }
    main=paste(colnames(L),main)
    if (length(main)==1) main=rep(main,ncol(L)) 
  if (x@probBin)
  {
     if (is.null(colRamp)) colRamp=.blueWhiteRed(21)
     for (k in 1:ncol(L)) 
     {
	hlp=L[,k]%*%getBasis(x)
        hlpCol=.binColsCont(hlp,colorLim,colRamp)
        if (is.null(set)) setL=which(L[,k]!=0) else setL=set
	.myFpDotPlot(x@fp,x@fmod,x@fset,hlpCol$binCol,hlpCol$plotOrder,set=setL,param=param,ask=ask,main=main[k],...)
        .myFpContPlot(c(hlp),hlpCol$binCol,type="l",main=main[k],ylab="Fingerprint-contrast",...)        	
     }
  
  } else
  {
    if (is.null(colRamp)) colRamp=.jet.colors2(512)
    colorLim=colorLim*c(-1,1)
    for (k in 1:ncol(L))
    .plotCont(t(L[,k]%*%getBasis(x)),nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,contour=contour,contourHlp=t((L[,k]%*%getBasis(x))),main=main[k],...)
  }
 }
})

#' @param disc selects the Principal Components or Discriminants that will be visualized
#' @param nRound optional to refine plots (axis labels of percentage of variability (discrimination) explained by PC (Discriminant)
#' @param plotType different plotTypes can be used. plotType="pcaPlot" plots the data in the principal component space using standard biplots, plotType="pcaInt" produces loading plots for interpretating principal components in the original space, plotType="pcaCont" produces plots to interpret contrasts in the PCA space, plotType="discPlot" produces discriminant biplots, plotType="discPlot2" produces stact discriminant plots similar to the De Roy et al. (2012) paper, plotType="discInt" produces loading plots for interpretating discriminants in the original space and plotType="discCont" produces  plots to interpret contrasts in the discriminant space with respect to features of the original bivariate projections of the multivariate distributions.
# @param fBasis flowBasis-class object used to construct flowPca object. Has to be supplied when using plotType="pcaInt", plotType="pcaCont", plotType="discInt" and plotType="discCont".
#' @param colRamp user defined colorRamp pallette can be used when producing interpertation plots for flowPca or flowDa objects 
#' @param colorLim optional argument for refining interpertation plots for flowPca or flowDa objects 
#' @param xAt optional to refine plots
#' @param gLab optional to refine plots
#' @param plotBox optional to refine plots
#' @param groups factor with group labels. Optional.
#' @examples
#' #
#' #
#' ###########################################
#' ###flowDa plots
#' ###########################################
#'
#' if(require(flowFDAExampleData)){
#' #Construct fPca object 
#' fPca=flowPca(fbasis)
#'
#' #Make plots
#' plot(fPca,disc=1:3)
#' plot(fPca,plotType="pcaInt",disc=1)
#' }
#' @rdname plot-methods
#' @aliases plot,flowPca-method

setMethod("plot","flowPca",function(x,y,disc=1:2,nRound=1,plotType="pcaPlot",fBasis=NULL,L=NULL,pow=.75,ask=TRUE,colRamp=NULL,colorLim=NULL,contour=FALSE,xAt=NULL,gLab=NULL,plotBox=TRUE,groups=as.factor(rep(1,nSet(x))),main=NULL,set=NULL,...)
{
  if (max(disc)>ncol(getPca(x)$x)) stop("Redefine disc argument. Less discriminants available than in the disc argument.")
  
  if (plotType=="pcaPlot")
  { if (length(disc)>=2) .plotScores(getPca(x)$x[,disc],paste("PC",disc," (",round(getPca(x)$sdev[disc]^2/sum(getPca(x)$sdev^2)*100,nRound),"%)",sep=""),col=as.double(groups),pch=as.double(groups),main=main,...)
    else .plotScores(cbind(getPca(x)$x[,disc],as.double(groups)),labels=c(paste("PC",disc," (",round(getPca(x)$sdev[disc[1]]^2/sum(getPca(x)$sdev^2)*100,nRound),"%)",sep=""),""),yaxt="n",col=as.double(groups),pch=as.double(groups),main=main,...)
  }
 
  if (plotType=="plotPca2")
{
.pcaPlot2(getPca(x),groups,xAt=xAt,gLab=gLab,plotBox=plotBox,disc=disc,main=main,...)
}

 
  if (plotType=="pcaInt")
  { 
   if (length(disc)>1) { warning("more than 1 discriminant is given, only the first one is used")
   disc=1}
   if (x@probBin){
    if(is.null(fBasis)) stop("fBasis object is not defined")
     if (is.null(colRamp)) colRamp=.blueWhiteRed(51)
     rot=getPca(x)$rotation[,disc]
     hlp<-.binColsCont(rot,max(abs(rot)),colRamp)
     par(ask=ask)
     param=getParam(x)
     for (i in 1:nrow(param))
     plot(fBasis@fmod,bin_col=hlp$binCol[order(hlp$plotOrder)],parameters=c(which(colnames(fBasis@fset)==param[i,1]),which(colnames(fBasis@fset)==param[i,2])),showbins=order(hlp$plotOrder),main=main,...)
     .myFpContPlot(rot,hlp$binCol,type="l",main=main,ylab="PC Loading",...)
   } else
   {if (is.null(colRamp)) colRamp=.jet.colors2(512)
    .plotCont(getPca(x)$rotation[,disc],nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,main=main,...)}
  }
  
   if (plotType=="pcaCont")
  {
 if (length(disc)>1) {warning("more than 1 discriminant is given, only the first one is used")
  disc=disc[1]}
if(class(L) !="matrix") L<-matrix(L)
    if(is.null(fBasis)) stop("fBasis object is not defined")
    if (nrow(L)!=nrow(getBasis(fBasis))) stop("contrast L has the wrong dimension")
    if (is.null(colorLim))
	{ for (k in 1:ncol(L))
          colorLim=max(abs(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,disc]))),colorLim)
        }
    main=paste(colnames(L),main)
    if (length(main)==1) main=rep(main,ncol(L)) 

    if (fBasis@probBin)
    {
    	if (is.null(colRamp)) colRamp=.blueWhiteRed(21)
    	for (k in 1:ncol(L))
	{         
         rot=t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,disc]))
        rotCol=.binColsCont(rot,colorLim)
        if (is.null(set)) setL=which(L[,k]!=0) else setL=set
	.myFpDotPlot(fBasis@fp,fBasis@fmod,fBasis@fset,rotCol$binCol,rotCol$plotOrder,set=setL,param=param,ask=ask,main=main,...)
	.myFpContPlot(c(L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis)))),rotCol$binCol,type="l",main=main,ylab="Fingerprint-contrast",...)
	}
    } else
    {
        if (is.null(colRamp)) colRamp=.jet.colors2(512)
	colorLim=colorLim*c(-1,1)
    	main=paste(colnames(L),main)
    	if (length(main)==1) main=rep(main,k) 
    	for (k in 1:ncol(L))
    	.plotCont(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,disc])),nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,contour=contour,contourHlp=t((L[,k]%*%getBasis(fBasis))),main=main[k],...)
    }
   }
})

#' @param cex.Axis is used to refine flowDa plots
#' @param cex.Text can be used to refine flowDa plots
#' @examples
#' #
#' #
#' ###########################################
#' ###flowDa plots
#' ###########################################
#' if(require(flowFDAExampleData)){
#' data(group)
#' data(fbasis)
#' #Construct flowDa object
#' fDa=flowDa(fbasis,groups=group,nPca=.95)
#' #Make plots, the flowDa class extends flowPca (PCA is performed before discriminant analysis). Hence all plots for flowPca objects can be made for flowDa objects. 
#' plot(fDa,disc=1:3,plotType="pcaPlot")
#' plot(fDa,disc=1:3)
#' plot(fDa,disc=1:3,plotType="discPlot2")
#' plot(fDa,disc=1,plotType="discInt")
#' }
#' @rdname plot-methods
#' @aliases plot,flowDa-method
setMethod("plot","flowDa",function(x,y,disc=1:2,nRound=1,plotType="discPlot",fBasis=NULL,L=NULL,pow=.75,ask=TRUE,colRamp=NULL,colorLim=NULL,contour=FALSE,xAt=NULL,gLab=NULL,plotBox=TRUE,groups=as.factor(rep(1,nSet(x))),cex.Axis=.7,cex.Text=1,cex=1,main=NULL,set=NULL,...)
{
  if (max(disc)>ncol(getDa(x)$M)) stop("Redefine disc argument. Less discriminants available than in the disc argument.")
  
  if (plotType=="discPlot")
  { if (length(disc)>=2) .plotScores(getDa(x)$M[,disc],paste("D",disc," (",round(getDa(x)$lambda[disc]/sum(getDa(x)$lambda)*100,nRound),"%)",sep=""),col=as.double(getGroups(x)),pch=as.double(getGroups(x)),main=main,cex=cex,...)
    else .plotScores(cbind(getDa(x)$M[,disc],as.double(getGroups(x))),labels=c(paste("Discriminant ",disc," (",round(getDa(x)$lambda[disc[1]]/sum(getDa(x)$lambda)*100,nRound),"%)",sep=""),""),yaxt="n",col=as.double(getGroups(x)),pch=as.double(getGroups(x)),main=main,cex=cex,...)
  }
if (plotType=="discPlot2")
{.fdaPlot2(fdaObj=getDa(x),exp=getGroups(x),disc=disc,gLab=gLab,cex.Axis=cex.Axis,cex.Text=cex.Text,xAt=xAt,plotBox=plotBox,nRound=nRound,main=main,cex=cex,...)
}
  if (plotType=="pcaPlot")
  { if (length(disc)>=2) .plotScores(getPca(x)$x[,disc],paste("PC",disc," (",round(getPca(x)$sdev[disc]^2/sum(getPca(x)$sdev^2)*100,nRound),"%)",sep=""),col=as.double(groups),pch=as.double(groups),main=main,cex=cex,...)
    else .plotScores(cbind(getPca(x)$x[,disc],as.double(groups)),labels=c(paste("PC",disc," (",round(getPca(x)$sdev[disc[1]]^2/sum(getPca(x)$sdev^2)*100,nRound),"%)",sep=""),""),yaxt="n",col=as.double(groups),pch=as.double(groups),main=main,cex=cex,...)
  }
 
  if (plotType=="plotPca2")
{
.pcaPlot2(getPca(x),groups,xAt=xAt,gLab=gLab,plotBox=plotBox,disc=disc,main=main,cex=cex,...)
}

 
 
  if (plotType=="pcaInt")
  { 
   if (length(disc)>1) { warning("more than 1 discriminant is given, only the first one is used")
   disc=1}
   if (x@probBin) 
   { if(is.null(fBasis)) stop("fBasis object is not defined")
     if (is.null(colRamp)) colRamp=.blueWhiteRed(51)
     rot=getPca(x)$rotation[,disc]
     hlp<-.binColsCont(rot,max(abs(rot)),colRamp)
     par(ask=ask)
     param=getParam(x)
     for (i in 1:nrow(param))
plot(fBasis@fmod,bin_col=hlp$binCol[order(hlp$plotOrder)],parameters=c(which(colnames(fBasis@fset)==param[i,1]),which(colnames(fBasis@fset)==param[i,2])),showbins=order(hlp$plotOrder),main=main,...)
     .myFpContPlot(rot,hlp$binCol,type="l",main=main,...)
    } else
    {
    if (is.null(colRamp)) colRamp=.jet.colors2(512)
    .plotCont(getPca(x)$rotation[,disc],nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,main=main,...)
    }
  }
  
   if (plotType=="pcaCont")
  {
 if (length(disc)>1) {warning("more than 1 discriminant is given, only the first one is used")
  disc=disc[1]}
if(class(L) !="matrix") L<-matrix(L)
    if(is.null(fBasis)) stop("fBasis object is not defined")
    if (nrow(L)!=nrow(getBasis(fBasis))) stop("contrast L has the wrong dimension")
    if (is.null(colorLim))
	{ for (k in 1:ncol(L))
          colorLim=max(abs(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,disc]))),colorLim)
        }
    main=paste(colnames(L),main)
    if (length(main)==1) main=rep(main,ncol(L)) 

    if (fBasis@probBin)
    {
    	if (is.null(colRamp)) colRamp=.blueWhiteRed(21)
    	for (k in 1:ncol(L))
	{         
         rot=t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,disc]))
        rotCol=.binColsCont(rot,colorLim)
        if (is.null(set)) setL=which(L[,k]!=0) else setL=set
plot(fBasis@fmod,bin_col=hlp$binCol[order(hlp$plotOrder)],parameters=c(which(colnames(fBasis@fset)==param[i,1]),which(colnames(fBasis@fset)==param[i,2])),showbins=order(hlp$plotOrder),main=main,...)
     .myFpContPlot(rot,hlp$binCol,type="l",main=main,ylab="Fingerprint-contrast",...)
	}
    } else
    {
        if (is.null(colRamp)) colRamp=.jet.colors2(512)
	colorLim=colorLim*c(-1,1)
    	main=paste(colnames(L),main)
    	if (length(main)==1) main=rep(main,k) 
    	for (k in 1:ncol(L))
    	.plotCont(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,disc])),nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,contour=contour,contourHlp=t((L[,k]%*%getBasis(fBasis))),main=main[k],...)
    }
   }


  if (plotType=="discInt")
  {
    if (length(disc)>1) {warning("more than 1 discriminant is given, only the first one is used")
    disc=disc[1]}
    if (x@probBin)
    {
    if(is.null(fBasis)) stop("fBasis object is not defined")
     if (is.null(colRamp)) colRamp=.blueWhiteRed(51)
     rot=getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc]
     hlp<-.binColsCont(rot,max(abs(rot)),colRamp)
     par(ask=ask)
     param=getParam(x)
     for (i in 1:nrow(param))
     plot(fBasis@fmod,bin_col=hlp$binCol[order(hlp$plotOrder)],parameters=param[i,],showbins=order(hlp$plotOrder),xlab=param[i,1],ylab=param[i,2],main=main,...)
     .myFpContPlot(rot,hlp$binCol,type="l",main=main,ylab="Discriminant loading",...)
    } else
   { if (is.null(colRamp)) colRamp=.jet.colors2(512)
    .plotCont(getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc],nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,main=main,cex=cex,...)
}
  }
  if (plotType=="discCont")
  {
 if (length(disc)>1) {warning("more than 1 discriminant is given, only the first one is used")
  disc=disc[1]}
if(class(L) !="matrix") L<-matrix(L)
    if(is.null(fBasis)) stop("fBasis object is not defined")
    if (nrow(L)!=nrow(getBasis(fBasis))) stop("contrast L has the wrong dimension")
    if (is.null(colorLim))
	{ for (k in 1:ncol(L))
          colorLim=max(abs(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc]))),colorLim)
        }
    main=paste(colnames(L),main)
    if (length(main)==1) main=rep(main,ncol(L)) 
    if (fBasis@probBin)
    {
    	if (is.null(colRamp)) colRamp=.blueWhiteRed(21)
 	for (k in 1:ncol(L))
	{         
         rot=t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc]))        
         rotCol=.binColsCont(rot,colorLim)
        if (is.null(set)) setL=which(L[,k]!=0) else setL=set
	.myFpDotPlot(fBasis@fp,fBasis@fmod,fBasis@fset,rotCol$binCol,rotCol$plotOrder,set=setL,param=param,ask=ask,main=main[k])
	.myFpContPlot(c(L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis)))),rotCol$binCol,type="l",main=main[k],ylab="Fingerprint-contrast",...)
	}
    } else
    {
        if (is.null(colRamp)) colRamp=.jet.colors2(512)
	colorLim=colorLim*c(-1,1)
    	for (k in 1:ncol(L))
    .plotCont(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc])),nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,contour=contour,contourHlp=t((L[,k]%*%getBasis(fBasis))),main=main[k],cex=cex,...)
    }
#    if(class(L) !="matrix") L<-matrix(L)
#    if(is.null(fBasis)) stop("fBasis object is not defined")
#    if (is.null(colRamp)) colRamp=.jet.colors2(512)
#    if (nrow(L)!=nrow(getBasis(fBasis))) stop("contrast L has the wrong dimension")
#    if (is.null(colorLim))
#	{ for (k in 1:ncol(L))
#          colorLim=max(abs(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc]))),colorLim)
#        }
#    colorLim=colorLim*c(-1,1)
#    main=paste(colnames(L),main)
#    if (length(main)==1) main=rep(main,k) 
#    for (k in 1:ncol(L))
#    .plotCont(t((L[,k]%*%(getBasis(fBasis)-matrix(1,nrow=nrow(getBasis(fBasis)),ncol=1)%*%colMeans(getBasis(fBasis))))*c(getPca(x)$rotation[,1:nPca(x)]%*%getDa(x)$delta[,disc])),nbin=nbin(x),param=getParam(x),ask=ask,colRamp=colRamp,pow=pow,colorLim=colorLim,contour=contour,contourHlp=t((L[,k]%*%getBasis(fBasis))),main=main[k],cex=cex,...)
#  
 }
})

