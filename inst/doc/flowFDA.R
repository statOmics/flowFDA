### R code from vignette source 'flowFDA.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: flowFDA.Rnw:99-108
###################################################
library(flowFDA)
#fset<-read.flowSet(path="~/Dropbox/LabMet/flowcytometry/stress_test_2/",
#transformation=FALSE)
#fset
#
##subset feet to reduce memory footprint
#param=c("SS Log","FL 1 Log","FL 3 Log")
#fset=fset[,param]
#fset


###################################################
### code chunk number 2: flowFDA.Rnw:112-117
###################################################
#mytrans<-function(x) x/2^16
#fset<-transform("FL 1 Log"=mytrans,"FL 3 Log"=mytrans,"SS Log"=mytrans)%on%fset
#rg <- rectangleGate(filterId="myRectGate", list("SS Log"=c(1/2^17, Inf),
#"FL 1 Log"=c(1/2^17, Inf),"FL 3 Log"=c(1/2^17,Inf)))
#fset<-Subset(fset,rg)


###################################################
### code chunk number 3: flowFDA.Rnw:120-121
###################################################
#logtrans<- function(x) log(x)


###################################################
### code chunk number 4: flowFDA.Rnw:125-131
###################################################
#construct experiment factor 
#files<-list.files(path="~/Dropbox/LabMet/flowcytometry/stress_test_2/",pattern=".fcs")
#expHlp<-unlist(strsplit(files,split="_replicate")) 
#dim(expHlp)<-c(2,length(fset))
#group<-as.factor(expHlp[1,])
#nGroup<-nlevels(group)


###################################################
### code chunk number 5: CompanionPkg (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite("flowFDAExampleData")


###################################################
### code chunk number 6: flowFDA.Rnw:145-147
###################################################
library(flowFDAExampleData)
library(flowFDA)


###################################################
### code chunk number 7: flowFDA.Rnw:150-156
###################################################
data(fset)
data(group)
param=c("SS Log","FL 1 Log","FL 3 Log")
nGroup=nlevels(group)
nSamp=length(fset)
groupLevels=levels(group)


###################################################
### code chunk number 8: flowFDA.Rnw:168-170
###################################################
fbasis=flowBasis(fset,param,nbin=128, bw=0.01)
fbasis


###################################################
### code chunk number 9: flowFDA.Rnw:180-182
###################################################
par(mfrow=c(2,2))
plot(fbasis,ask=FALSE,samples=3)


###################################################
### code chunk number 10: flowFDA.Rnw:192-195
###################################################
par(mfrow=c(2,3))
plot(fbasis,ask=FALSE,samples=group==groupLevels[1],main=groupLevels[1])
plot(fbasis,ask=FALSE,samples=group==groupLevels[4],main=groupLevels[4])


###################################################
### code chunk number 11: flowFDA.Rnw:206-212
###################################################
par(mfrow=c(2,2))
L=rep(0,length(group))
L[group==groupLevels[1]]=-1/sum(group==groupLevels[1])
L[group==groupLevels[4]]=1/sum(group==groupLevels[4])
plot(fbasis,L=L,ask=FALSE,main=paste(groupLevels[4],"-",groupLevels[1],sep=""),
contour=TRUE,contourLwd=4,contourLevel=c(-.04,.04))


###################################################
### code chunk number 12: flowFDA.Rnw:229-238
###################################################
#construct flowPca object
fPca=flowPca(fbasis) 

#perform model based clustering, 
#use n PCs so as to capture at least 95 % of the variability
nPca(fPca)<-.95
nPca(fPca) #number of PCs used for model based clustering
setClust(fPca)<-Mclust(getPcaScore(fPca,nPca(fPca))) #Model based clustering
cbind(as.character(getClustClass(fPca)),as.character(group)) # cluster class labels and real grouping


###################################################
### code chunk number 13: flowFDA.Rnw:242-246
###################################################
par(mfrow=c(1,2))
plot(fPca,groups=getClustClass(fPca),main="Kernel Dens. (Clustering)") 
plot(fPca,groups=group,main="Kernel Dens. (Treatment)")
legend("topleft",legend=c("c","h24","h3","n24","n3"),pch=1:5,col=1:5) 


###################################################
### code chunk number 14: flowFDA.Rnw:274-288
###################################################
intSamples=3 #for the group average of first group set intSamples=which(group=groupLevels[1])
par(mfrow=c(2,2))
plot(fPca,groups=group,main="Treatment")
pcX=mean(getPca(fPca)$x[intSamples,1])
pcY=mean(getPca(fPca)$x[intSamples,2])
arrows(x0=pcX,x1=pcX,y0=-2,y1=pcY)

 #PCA is done after centering
 # interpretation in terms of contrast to average bivariate density
 #contrast between average bivariate density of intSamples vs overall average
L=rep(-1/nSamp,nSamp)
L[intSamples]=L[intSamples]+1/length(intSamples) 

plot(fPca,fBasis=fbasis,disc=1,plotType="pcaCont",L=L,ask=FALSE,main="PC 1",contour=TRUE,contourLwd=3,contourLevel=c(-.04,.04))


###################################################
### code chunk number 15: flowFDA.Rnw:305-310
###################################################
#supervised, class labels are needed
#select first few PC's which explain more than 95% of the variability in the original fingerprint. 
#####Discriminant analysis for kernel dens.
fDa=flowDa(fbasis,groups= group, nPca=.95) 
fDa


###################################################
### code chunk number 16: flowFDA.Rnw:316-320
###################################################
par(mfrow=c(1,2))
plot(fDa,groups=group,main="Kernel Dens. PCA",plotType="pcaPlot")
plot(fDa,main="Kernel Dens. DA")
legend("bottomleft",legend=c("c","h24","h3","n24","n3"),pch=1:5,col=1:5)


###################################################
### code chunk number 17: flowFDA.Rnw:343-348
###################################################
nPerm=100 
#Only 100 permutations are used 
#so as to restrict the computational burden when generating the vignette
disc=1:2 #Test only in the space of first 2 discriminants
fDa=flowDaTest(fDa,disc=disc,nPerm)


###################################################
### code chunk number 18: flowFDA.Rnw:352-355
###################################################
data(fDa)
adjustedPvalues=pAdjustMx(getMpc(fDa)$pValuePerm)
adjustedPvalues


###################################################
### code chunk number 19: flowFDA.Rnw:376-386
###################################################
nSamp=nSet(fDa)
L<-rep(0,nSamp)
L[group==groupLevels[4]]<-1/sum(group==groupLevels[4])
L[group==groupLevels[1]]<--1/sum(group==groupLevels[1])
par(mfrow=c(2,2))
plot(fDa)
legend("bottomleft",legend=c("c","h24","h3","n24","n3"),pch=1:5,col=1:5)
disc=1
plot(fDa,fBasis=fbasis,L=L,ask=FALSE,plotType="discCont",disc=disc,
contour=TRUE,contourLevel=c(-.04,.04),contourLwd=4)


###################################################
### code chunk number 20: flowFDA.Rnw:397-407
###################################################
nSamp=nSet(fDa)
L<-rep(0,nSamp)
L[group==groupLevels[2]]<-1/sum(group==groupLevels[2])
L[group==groupLevels[1]]<--1/sum(group==groupLevels[1])
par(mfrow=c(2,2))
plot(fDa)
legend("bottomleft",legend=c("c","h24","h3","n24","n3"),pch=1:5,col=1:5)
disc=2
plot(fDa,fBasis=fbasis,L=L,ask=FALSE,plotType="discCont",disc=disc,
contour=TRUE,contourLevel=c(-.04,.04),contourLwd=4)


###################################################
### code chunk number 21: flowFDA.Rnw:418-502
###################################################
library(flowFDAExampleData)
library(flowFDA)
data(fbasis)
data(fDa)
nSamp=nSet(fDa)
group=getGroups(fDa)
nGroup=nlevels(group)
groupLevels=levels(group)
sampleNames=rownames(getBasis(fbasis))

##############################
#Generate original Fingerprint plots
##############################

#uncomment to plot all bivariate distributions
#par(mfrow=c(1,3))
#for (i in 1:nSamp) plot(fbasis,sample=i,ask=TRUE,main=sampleNames[i])

#uncomment to create a pdf with all bivariate distributions plots
#pdf("allBasis.pdf",height=7,width=15)
#par(mfrow=c(1,3))
##cex to enlarge font
#for (i in 1:nSamp) 
#plot(fbasis,sample=i,ask=FALSE,main=sampleNames[i],
#cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
#dev.off()


#uncomment to plot all average bivariate distributions for each group
#par(mfrow=c(1,3))
#for (i in 1:nGroup) plot(fbasis,sample=group==groupLevels[i],ask=TRUE,main=groupLevels[i])

#uncomment to create a pdf with all average bivariate distributions plots for each group
#pdf("allGroupBasis.pdf",height=7,width=15)
#par(mfrow=c(1,3))
#for (i in 1:nGroup) 
#plot(fbasis,sample=group==groupLevels[i],ask=FALSE,main=groupLevels[i],
#cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
#dev.off()

##############################################
#Contrast interpretation plots
##############################################

#extract groups involved in contrasts from p value object
comp=strsplit(rownames(getMpc(fDa)$p),split="-") 

#create all corresponding contrasts
L=sapply(comp,function(x,group) 
(group==x[1])/sum(group==x[1])-(group==x[2])/sum(group==x[2]),group=group)
colnames(L)<-rownames(getMpc(fDa)$p)

#uncomment to generate contrast plot of flowBasis Object
#par(mfrow=c(1,3))
#plot(fbasis,L=L,ask=TRUE,contour=TRUE,contourLwd=4,contourLevel=c(-.04,.04))

#uncomment to generate a pdf of the contrast plot in original space
#pdf("contrastInterpretationPlotsBasis.pdf",height=7,width=15)
#par(mfrow=c(1,3))
#plot(fbasis,L=L,ask=FALSE,contour=TRUE,contourLwd=4,contourLevel=c(-.04,.04),
#cex.main=1.5,cex.axis=1.5,cex.lab=1.5)
#dev.off()

adjustedPvalues=pAdjustMx(getMpc(fDa)$p)
#uncomment to make discriminant interpretation plots in plot window
#par(mfrow=c(1,3))
#disc=1
#plot(fDa,fBasis=fbasis,L=L,ask=TRUE,plotType="discCont",disc=disc,contour=TRUE,
#contourLevel=c(-.04,.04),contourLwd=4,
#main=paste("\n D", disc," p=",round(adjustedPvalues[,disc],3),sep="")) 
#disc=2
#plot(fDa,fBasis=fbasis,L=L,ask=TRUE,plotType="discCont",disc=disc,contour=TRUE,
#contourLevel=c(-.04,.04),contourLwd=4,
#main=paste("\n D", disc," p=",round(adjustedPvalues[,disc],3),sep="")) 

#uncomment to create a pdf of the discriminant interpretation plots of contrasts
#pdf("contrastInterpretationPlotsDiscriminant.pdf",height=7,width=15)
#par(mfrow=c(1,3))
#for (disc in 1:2)
#plot(fDa,fBasis=fbasis,L=L,ask=FALSE,plotType="discCont",disc=disc,contour=TRUE,
#contourLevel=c(-.04,.04),contourLwd=4,
#main=paste("\n D", disc," p=",round(getMpc(fDa)$p[,disc],3),sep=""),
#cex.main=1.5,cex.axis=1.5,cex.lab=1.5) 
#dev.off()


