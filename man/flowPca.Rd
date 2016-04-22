\name{flowPca}
\alias{flowPca}
\title{flowPca constructor}
\usage{
  flowPca(fbasis)
}
\arguments{
  \item{fbasis}{A \code{\link{flowBasis-class}} object for
  which a discriminant basis is desired.}
}
\value{
  An instance of an object of the
  \code{\link{flowPca-class}}
}
\description{
  Constructor for \code{\link{flowPca-class}}
}
\details{
  A flowPca object is a reduced representation of a
  \code{\link{flowBasis-class}} and a
  flowCore:flowSet-class with respect to the bivariate
  distributions of the flow channels that are considered.
  Consists of a principal component decomposition of the
  coefficients of the estimated bivariate kernel densities.
  The PCs can be interpreted with respect to the features
  in the original bivariate distributions.
}
\examples{
# load a flowSet to use as an example basis object
if(require(flowFDAExampleData)){
data(fbasis)
data(group)
data(param)

flowpca<-flowPca(fbasis)
plot(flowpca)

####Plot of the loadings the first PC of the daobject
pc=1
plot(flowpca,disc=pc,main=paste("PC",pc),plotType="pcaInt")

#interpretation of score on PC1 for first flowcytometry sample
#pc=1
#L=rep(0,length(group))
#L[1]<-1
#plot(flowpca,flowbasis,disc=pc,L=L,plotType="pcaCont")
}
}
\references{
  De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N.
  (2012). Flow cytometry for fast microbial community
  fingerprinting. Water Research, 46 (3), 907-919.
}
\seealso{
  \code{\link{flowBasis-class}},
  \code{\link{flowPca-class}} and
  \code{\link{flowDa-class}}
}

