###ColorRamps
.jet.colors <-colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
.jet.colors2 <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","white", "yellow", "#FF7F00", "red", "#7F0000"))
.blueBlackRed<-colorRampPalette(c("blue","black","red"))
.blueWhiteRed<-colorRampPalette(c("darkblue","blue","white","red","darkred"))
.blueGreyRed<-colorRampPalette(c("blue","#BFBFBF","red"))

.blueWhiteRed<-colorRampPalette(c("blue","white","red"))
.greenWhiteOrange<-colorRampPalette(c("green","white","orange"))

#########print methods
print.Mclust<-function(x, digits = getOption("digits"), ...)
{
if (as.character(x$call)[1]=="<undef>") cat("cluster object has not been initalized yet\n") else mclust::print.Mclust(x, digits = digits,...)
}

print.flowBasis<-function(x,...)
{
cat("flowBasis object\n")
  if (x@probBin) cat("Probability Binning using",nbin(x),"bins.\n\n") else cat("Kernel Density Estimation on grid of",nbin(x),"x",nbin(x),"\nKernel density bandwith:",getBw(x),"\n\n")        
  cat("bivariate densities for channels\n")
  print(getParam(x))
}
print.flowDa<-function(x,digits=5,...)
{
  cat("flowDa object\n")
if (x@probBin) cat("Probability Binning using",nbin(x),"bins.\n\n") else cat("Kernel Density Estimation on grid of",nbin(x),"x",nbin(x),"\n\n")        
  cat("channels\n")
  print(getParam(x))
  if (length(getMpc(x))>0)
  {
    cat("\ntests\n")
    mpcPrint=cbind(getMpc(x)$contrasts,getMpc(x)$pValuePerm)
    colnames(mpcPrint)<-paste(c(rep("",ncol(mpcPrint)/2),rep("p-value ",ncol(mpcPrint)/2)),colnames(mpcPrint),sep="")
    print(mpcPrint,digits=digits)
  }
}

print.flowPca<-function(x,...)
{
  cat("flowPca object\n")
if (x@probBin) cat("Probability Binning using",nbin(x),"bins.\n\n") else cat("Kernel Density Estimation grid of",nbin(x),"x",nbin(x),"\n\n")        
  cat("channels\n")
  print(getParam(x))
}

setMethod("show","flowBasis",function(object) print.flowBasis(object))
setMethod("show","flowDa",function(object) print.flowDa(object))
setMethod("show","flowPca",function(object) print.flowPca(object))

.mpcAll<-function(M,groups,mpc,disc,what) {
if(!is.element(what,c("statistic","p.value"))) stop("what-argument should be statistic or p.value")
tests<-t(sapply(1:nrow(mpc),function(x,M,mpc,groups,disc,what) apply(M[,disc],2,function(x,ind,mpc,groups) t.test(x[groups==mpc[ind,1]],x[groups==mpc[ind,2]])[[what]],mpc=mpc,groups=groups,ind=x),M=M,mpc=mpc,groups=groups,disc=disc,what=what))
#colnames(tests)<-paste("D",disc,sep="")
rownames(tests)<-paste(mpc[,2],mpc[,1],sep="-")
return(tests)
}

.flowDaCvHlp<-function(fbasis,indVal,groups,bw,nPca=length(indVal),disc=1:2,nbin=128,normalize=function(x) x/max(x))
{
	ndisc<-length(disc)
	fBasisCal<-fBasisVal<-fbasis
	fBasisCal@basis<-getBasis(fbasis)[-indVal,]
	fBasisVal@basis<-getBasis(fbasis)[indVal,]	
	fDaCal<-flowDa(fBasisCal,groups[-indVal],nPca=nPca)
	groupMeans<-lm(getDa(fDaCal)$M~-1+getGroups(fDaCal))$coef
	groupLevels=levels(groups)
	nGroup<-nlevels(groups)
	return(sapply(1:length(indVal),function(k) (1:nGroup)[order(apply(groupMeans,1,function(x,test,disc,ndisc) sum((test[disc[1:ndisc]]-x[disc[1:ndisc]])^2),test=(getBasis(fBasisVal)[k,]-getPca(fDaCal)$center)%*%getPca(fDaCal)$rotation[,1:nPca]%*%getDa(fDaCal)$delta,ndisc=ndisc,disc=disc))[1]]))
}

setMethod("flowDaCv",signature="flowDa",definition=function(fset,param,cv,groups,nbins=128,pcs=2:15,bw=0.05,file=NULL)
{
	nCv<-nrow(cv)
	fbasis<-flowBasis(fset,param,nbins,bw)
	evalCv<-lapply(pcs,function(xx,fbasis,cv,groups,file) {evalInner<-sapply(1:nCv,function(x,fbasis,cv,groups,bw,nPca) {
	cat("bw",bw,"pca",nPca,"cv",x,"/",nrow(cv),"\n")
	.flowDaCvHlp(fbasis,cv[x,],groups,bw,nPca)
	},fbasis=fbasis,cv=cv,groups=groups,bw=bw,nPca=xx)
	if (!is.null(file)) cat(evalInner,file=paste(file,xx,".txt",sep=""))
	return(evalInner)},fbasis=fbasis,cv=cv,groups=groups,file=file)
	#names(evalCv)<-paste("npc",pcs,sep="")	
	return(list(evalSet=evalCv,cvSet=sapply(1:nCv,function(x,cv,groups) as.double(groups)[cv[x,]],cv=cv,groups=groups)))
})



#####plot methods
#Low Level plot functions
.plotScores<-function(scores,labels=NULL,...)      
{
  scores<-as.matrix(scores)
  if (ncol(scores)>2) pairs(scores,labels,...) else plot(scores,xlab=labels[1],ylab=labels[2],...)
}

.plotCont<-function(hlp,nbin,param,pow=.75,ask=TRUE,colorLim=NULL,colRamp=NULL,contour=FALSE,contourHlp=NULL,contourCol=c("blue","red"),contourLevel=c(-0.02,0.02),contourLab=c("-","+"),contourLwd=2,...)
{
  hlp<-as.matrix(hlp)
  par(ask=ask)
  if(is.null(colRamp)) colRamp<-.jet.colors2(512)
  if (ncol(hlp)>1) 
  {
    warning("more than one component selected, compontents are averaged")
    hlp<-rowSums(hlp)
  }
  for (i in 1:(length(hlp)/nbin^2))
  {
    hlp2<-matrix(hlp[1:(nbin*nbin)+(i-1)*(nbin*nbin)],ncol=nbin)
    if (is.null(colorLim)) colorLim=max(abs(hlp))*c(-1,1)
    hlp2[nbin,nbin]<-colorLim[1]
    hlp2[nbin,nbin-1]<-colorLim[2]
    image(sign(hlp2)*abs(hlp2)^pow,xlab=param[i,1],ylab=param[i,2],col=colRamp,...)
    if (contour) 
    {
      hlp3<-matrix(contourHlp[1:(nbin*nbin)+(i-1)*(nbin*nbin)],ncol=nbin)
      contour(matrix(hlp3,ncol=nbin),add=TRUE,col=contourCol[1],level=contourLevel[sign(contourLevel)==-1],lwd=contourLwd,labels=contourLab[1])
      contour(matrix(hlp3,ncol=nbin),add=TRUE,col=contourCol[2],level=contourLevel[sign(contourLevel)==1],lwd=contourLwd,labels=contourLab[2])
    }
  }
}

.fdaPlot2<-function(fdaObj,exp,disc=1:3,gLab=NULL,cex.Axis=.7,cex.Text=1,xlim=NULL,xAt=NULL,plotBox=TRUE,nRound=1,...)
{
###disc is standardly equals 1:3 to plot first discriminants
###use disc=c(1,4) to plot discriminant 1 and 4 or disc=4 to only print discriminant 4
###change cex.axis if you want to alter font of axis
###change cex.text if you want to alter the fontsize of the labels
expGroup=levels(exp)
nG<-nlevels(exp)
nRep<-sapply(expGroup,function(x) sum(exp==x))
ndisc<-length(disc)
if (is.null(xlim))
{ 
	xlim=range(fdaObj$M)
	xlim[1]=xlim[1]-0.1*(xlim[2]-xlim[1])
}
plot(range(fdaObj$M),c(0,0),col=0,ylim=c(nG*ndisc+3*(ndisc-1)+2,1),axes=FALSE,ylab="",xlab="Score",xlim=xlim,...)
for (j in 1:ndisc) for (i in 1:nG) points(fdaObj$M[exp%in%expGroup[i],disc[j]],rep(i+(j-1)*nG+(j-1)*3+2,nRep[i]),col=i,pch=i)
axis(1,at=xAt)
gPos<-rep(0,nG*ndisc)
for (j in 1:ndisc) gPos[(1:nG)+(j-1)*nG]<-(1:nG)+(j-1)*nG+(j-1)*3+2
xpos=xlim[1]
if (!is.null(gLab)) {axis(2,at=gPos,labels=rep(gLab,ndisc),cex.axis=cex.Axis)
   } else text(rep(xpos,nG*ndisc),gPos,labels=expGroup,pos=4,col=rep(1:nG,ndisc),cex=cex.Text)

if (ndisc>1) abline(h=((1:(ndisc-1))*nG+3*(1:(ndisc-1))))
xpos<-mean(xlim)
if (ndisc>1) text(rep(xpos,ndisc),c((0:(ndisc-1))*nG+3*(1:(ndisc))-2),paste(rep("Discriminant ",ndisc),disc,rep(" (",ndisc),round(fdaObj$lambda[disc]/sum(fdaObj$lambda)*100,nRound),rep("%)",ndisc),sep="")) else text(xpos,1,paste(rep("Discriminant ",disc),1:ndisc,rep(" (",ndisc),round(fdaObj$lambda[disc]/sum(fdaObj$lambda)*100,nRound),rep("%)",ndisc),sep="")) 
if (plotBox) box()
}
.pcaPlot2<-function(pcaObj,exp,disc=1:3,gLab=NULL,cex.Axis=.7,cex.Text=1,xlim=NULL,xAt=NULL,plotBox=TRUE,nRound=1,...)
{
###disc is standardly equals 1:3 to plot first discriminants
###use disc=c(1,4) to plot discriminant 1 and 4 or disc=4 to only print discriminant 4
###change cex.axis if you want to alter font of axis
###change cex.text if you want to alter the fontsize of the labels

ndisc<-length(disc)
expGroup=levels(exp)
nG<-nlevels(exp)
nRep<-sapply(expGroup,function(x) sum(exp==x))
if (is.null(xlim))
{ 
	xlim=range(pcaObj$x)
	xlim[1]=xlim[1]-0.1*(xlim[2]-xlim[1])
}
plot(range(pcaObj$x),c(0,0),col=0,ylim=c(nG*ndisc+3*(ndisc-1)+2,1),axes=FALSE,ylab="",xlab="Score",xlim=xlim,...)
for (j in 1:ndisc) for (i in 1:nG) points(pcaObj$x[exp%in%expGroup[i],disc[j]],rep(i+(j-1)*nG+(j-1)*3+2,nRep[i]),col=i,pch=i)
axis(1,at=xAt)
gPos<-rep(0,nG*ndisc)
for (j in 1:ndisc) gPos[(1:nG)+(j-1)*nG]<-(1:nG)+(j-1)*nG+(j-1)*3+2
xpos=xlim[1]
if (!is.null(gLab)) {axis(2,at=gPos,labels=rep(gLab,ndisc),cex.axis=cex.Axis)
   } else text(rep(xpos,nG*ndisc),gPos,labels=expGroup,pos=4,col=rep(1:nG,ndisc),cex=cex.Text)

if (ndisc>1) abline(h=((1:(ndisc-1))*nG+3*(1:(ndisc-1))))
xpos<-mean(xlim)
if (ndisc>1) text(rep(xpos,ndisc),c((0:(ndisc-1))*nG+3*(1:(ndisc))-2),paste(rep("PC ",ndisc),disc,rep(" (",ndisc),round(pcaObj$sdev[disc]^2/sum(pcaObj$sdev^2)*100,nRound),rep("%)",ndisc),sep="")) else text(xpos,1,paste(rep("PC ",disc),1:ndisc,rep(" (",ndisc),round(pcaObj$sdev[disc]^2/sum(pcaObj$sdev^2)*100,nRound),rep("%)",ndisc),sep="")) 
if (plotBox) box()
}

#' Function that extends p.adjust to matrices of p-values
#'
#' @param pvals matrix with p-values
#' @param method method to correct for multiple testing, c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"). Details can be found in stats::p.adjust.
#' @return matrix with adjusted pvalues
#' 
#' @examples
#' if(require(flowFDAExampleData)){ 
#' data(fDa)
#' pAdjustMx(getMpc(fDa)$p)
#' }
#' @importFrom stats 
#' @export
#' @rdname pAdjustMx
pAdjustMx<-function(pvals,method="holm")
{
padj<-p.adjust(pvals,method="holm")
dim(padj)<-dim(pvals)
rownames(padj)<-rownames(pvals)
colnames(padj)<-colnames(pvals)
return(padj)
}





#########################################################
#old functions for compatibility with flowFP implementation is De Roy et al. 2012
#########################################################
.binColsCont<-function(rotation,maxCol,colRamp=.blueWhiteRed(21),plotOrderAbs=TRUE)
{
	#Contrast: expT - expRef	
	#expT: factor level of treatment T
	#expRef: factor level of ref treatment
	#disc: discriminant for which color coding is calculated
        nPosCols=length(colRamp)%/%2
        binColHlp<-sign(rotation)*(abs(round(rotation/maxCol*nPosCols)))+nPosCols+1
        if (plotOrderAbs) return(list(binColHlp=binColHlp,plotOrder=abs(binColHlp-nPosCols-1),binCol=colRamp[binColHlp],nPosCols=nPosCols)) else return(list(binColHlp=binColHlp,plotOrder=binColHlp,binCol=colRamp[binColHlp],nPosCols=nPosCols))
}

testNames<-function(testAllObj)
{
out<-unlist(strsplit(names(testAllObj),split="\\."))
dim(out)<-c(4,length(out)/4)
return(t(out))
}

.myFpDotPlot<-function(fp,fmod,fset,binCols,plotOrder,set,param,nPosCol=max(plotOrder),ask=TRUE,xlim=NULL,ylim=NULL,...)
{
#make interpretation plot
#fp, flowFP object
#fmod flowFPModel
#fset flowSet
#binCols colors that represent the for the different bins
#plotOrder plot order for the bins, important for plots involving models with more than two variables as bins can be projected onto each other in 2D visualisation
#set index with data
#param, variables to plot
#nPosCol: number of positive colors that were used to construct the binCols color scheme

dotcols<-list()
dotOrder<-list()
tags = tags(fp)
for (i in 1:length(set)) dotcols[[i]]<-binCols[tags[[set[i]]]]
for (i in 1:length(set)) dotOrder[[i]]<-plotOrder[tags[[set[i]]]]
dotcols<-unlist(dotcols)
dotOrder<-unlist(dotOrder)
nPar=length(param)
par(ask=ask)
for (j in 1:(nPar-1))
for (k in (j+1):(nPar))
{
plothlp<-NULL
if (is.null(xlim)) xlim=range(unlist(fsApply(fset[set],function(x) range(x[,param[j]]))))
if (is.null(ylim)) ylim=range(unlist(fsApply(fset[set],function(x) range(x[,param[k]]))))
for (i in set)
if (is.null(plothlp)) plothlp<-exprs(fset[[i]])[,param[c(j,k)]] else  plothlp<-rbind(plothlp,exprs(fset[[i]])[,param[c(j,k)]])
plot(plothlp[dotOrder==0,],col=dotcols[dotOrder==0],pch=".",xlim=xlim,ylim=ylim,...)
for (i in 1:nPosCol) points(plothlp[dotOrder==i,],col=dotcols[dotOrder==i],pch=".")
}
}


.myFpContPlot<-function(Cont,binCols,...)
{
plot(Cont,...)
usr = par("usr")
color_wedge_h = (usr[4] - usr[3]) * 0.05
ypos = seq(from=usr[3] - color_wedge_h, to=usr[3] + color_wedge_h, length.out=2)
xpos = seq(from=0.5, to=length(binCols) + .5, length.out=length(binCols) + 1)
wedge = matrix(rep(c(1:length(binCols)), 1), ncol=length(binCols), nrow=1, byrow=TRUE)
image(xpos, ypos, t(wedge), col=binCols, add=TRUE, xlab=NA, ylab=NA, yaxt="n" )
abline(h=0,lty=2,col="grey")
}


