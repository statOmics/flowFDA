### R code from vignette source 'flowFDAProbabilityBinning.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: CompanionPkg (eval = FALSE)
###################################################
## source("http://www.bioconductor.org/biocLite.R")
## biocLite("flowFDAExampleData")


###################################################
### code chunk number 2: flowFDAProbabilityBinning.Rnw:106-114
###################################################
library(flowFDAExampleData)
library(flowFDA)
data(fset)
data(group)
param=c("SS Log","FL 1 Log","FL 3 Log")
nGroup=nlevels(group)
nSamp=length(fset)
groupLevels=levels(group)


###################################################
### code chunk number 3: flowFDAProbabilityBinning.Rnw:134-136
###################################################
fbasisPb=flowBasis(fset,param,nbin=128,probBin=TRUE)
fbasisPb


###################################################
### code chunk number 4: flowFDAProbabilityBinning.Rnw:146-148
###################################################
par(mfrow=c(2,2))
plot(fbasisPb,ask=FALSE,samples=3)


###################################################
### code chunk number 5: flowFDAProbabilityBinning.Rnw:158-161
###################################################
par(mfrow=c(2,4))
plot(fbasisPb,ask=FALSE,samples=group==groupLevels[1],main=groupLevels[1])
plot(fbasisPb,ask=FALSE,samples=group==groupLevels[4],main=groupLevels[4])


###################################################
### code chunk number 6: flowFDAProbabilityBinning.Rnw:174-180
###################################################
par(mfrow=c(2,2))
L=rep(0,length(group))
L[group==groupLevels[1]]=-1/sum(group==groupLevels[1])
L[group==groupLevels[4]]=1/sum(group==groupLevels[4])
par(mfrow=c(2,2))
plot(fbasisPb,L=L,set=which(L!=0),ask=FALSE,main=paste(groupLevels[4],"-",groupLevels[1],sep=""))


###################################################
### code chunk number 7: flowFDAProbabilityBinning.Rnw:189-198
###################################################
#construct flowPca object with probability binning basis
fPcaPb=flowPca(fbasisPb) 

#perform model based clustering, 
#use n PCs so as to capture at least 95 % of the variability
nPca(fPcaPb)<-.95
nPca(fPcaPb) #number of PCs used for model based clustering
setClust(fPcaPb)<-Mclust(getPcaScore(fPcaPb,nPca(fPcaPb))) #Model based clustering
cbind(as.character(getClustClass(fPcaPb)),as.character(group)) # cluster class labels and real grouping


###################################################
### code chunk number 8: flowFDAProbabilityBinning.Rnw:202-206
###################################################
par(mfrow=c(1,2))
plot(fPcaPb,groups=getClustClass(fPcaPb),main="Prob. Bin. (Clustering)") 
plot(fPcaPb,groups=group,main="Prob. Bin. (Treatment)")
legend("topleft",legend=c("c","h24","h3","n24","n3"),pch=1:5,col=1:5) 


###################################################
### code chunk number 9: flowFDAProbabilityBinning.Rnw:220-235
###################################################
intSamples=3 #for the group average of first group set intSamples=which(group=groupLevels[1])
layout(matrix(c(1,2,3,1,4,5),nrow=2,byrow=TRUE))
par(pty="s")
plot(fPcaPb,groups=group,main="Actual grouping")
pcX=mean(getPca(fPcaPb)$x[intSamples,1])
pcY=mean(getPca(fPcaPb)$x[intSamples,2])
arrows(x0=pcX,x1=pcX,y0=-4,y1=pcY)
intSamples=3 #for the group average of first group set intSamples=which(group=groupLevels[1])

 #PCA is done after centering
 # interpretation in terms of contrast to average bivariate density
 #contrast between average bivariate density of intSamples vs overall average
L=rep(-1/nSamp,nSamp)
L[intSamples]=L[intSamples]+1/length(intSamples) 
plot(fPcaPb,fBasis=fbasisPb,disc=1,plotType="pcaCont",L=L,ask=FALSE,main="PC 1")


###################################################
### code chunk number 10: flowFDAProbabilityBinning.Rnw:244-247
###################################################
#####Discriminant analysis for prob. bin.
fDaPb=flowDa(fbasisPb,groups= group, nPca=.95) 
fDaPb


###################################################
### code chunk number 11: flowFDAProbabilityBinning.Rnw:252-256
###################################################
par(mfrow=c(1,2))
plot(fDaPb,groups=group,main="Prob. Bin. PCA",plotType="pcaPlot")
plot(fDaPb,main="Prob. Bin. DA")
legend("topleft",legend=c("c","h24","h3","n24","n3"),pch=1:5,col=1:5)


###################################################
### code chunk number 12: flowFDAProbabilityBinning.Rnw:269-275
###################################################
nPerm=100 
#Only 100 permutations are used 
#so as to restrict the computational burden when generating the vignette
#nPerm=10000
disc=1:2 #Test only in the space of first 2 discriminants
fDaPb=flowDaTest(fDaPb,disc=disc,nPerm)


###################################################
### code chunk number 13: flowFDAProbabilityBinning.Rnw:280-282
###################################################
adjustedPvalues=pAdjustMx(getMpc(fDaPb)$pValuePerm)
adjustedPvalues


###################################################
### code chunk number 14: flowFDAProbabilityBinning.Rnw:290-300
###################################################
groupLevels=levels(group)
nSamp=nSet(fDaPb)
L<-rep(0,nSamp)
L[group==groupLevels[4]]<-1/sum(group==groupLevels[4])
L[group==groupLevels[1]]<--1/sum(group==groupLevels[1])
layout(matrix(c(1,2,3,1,4,5),nrow=2,byrow=TRUE))
par(pty="s")
plot(fDaPb)
disc=1
plot(fDaPb,fBasis=fbasisPb,L=L,ask=FALSE,plotType="discCont",disc=disc)


