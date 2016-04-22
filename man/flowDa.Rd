\name{flowDa}
\alias{flowDa}
\title{flowDa constructor}
\usage{
  flowDa(fbasis, groups, nPca = 2)
}
\arguments{
  \item{fbasis}{A \code{\link{flowBasis-class}} object for
  which a discriminant basis is desired.}

  \item{groups}{A vector of the numeric or factor type with
  the class labels}

  \item{nPca}{The number of PCs or the fraction of the
  total variation to be used in the PCA dimension reduction
  prior to the discrimination analysis}
}
\value{
  An instance of an object of the
  \code{\link{flowDa-class}}
}
\description{
  Constructor for \code{\link{flowDa-class}}
}
\details{
  A flowDa object is a reduced representation of a
  \code{\link{flowBasis-class}} and a
  flowCore:flowSet-class with respect to the bivariate
  distributions of the flow channels that are considered.
  The discriminants are derived by adopting Fisher's Method
  on the first \code{nPca} principal components that were
  derived from the coefficients of the bivariate kernel
  density estimation. The discriminants can be interpreted
  with respect to the features in the original bivariate
  distributions.
}
\examples{
if(require(flowFDAExampleData)){
# load a flowSet to use as an example basis object load(fbasis)
data(fbasis)
data(group)
data(param)
flowda<-flowDa(fbasis,group,nPca=6)
plot(flowda)

####Plot of the loadings the first discriminant of the daobject
#discriminant=1
#plot(flowda,disc=discriminant,plotType="discInt",main=paste("Discriminant",discriminant),cex.axis=1.5,cex.lab=2,cex.main=2)

####interpretation plot using a vector L to assess the contribution of centroid of the group 1 on the first discrimant
####Other definitions of L support the interpretation of the score for an individual sample or for contrasts between groups etc.
#discriminant=1
#L=group==levels(group)[1]
#plot(flowda,fbasis,disc=discriminant,plotType="discCont",L=L,main=paste("Discriminant",discriminant),cex.axis=1.5,cex.lab=2,cex.main=2)
}
}
\references{
  De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N.
  (2012). Flow cytometry for fast microbial community
  fingerprinting. Water Research, 46 (3), 907-919.
}
\seealso{
  \code{\link{flowBasis-class}},
  \code{\link{flowPca-class}}, \code{\link{flowDa-class}},
  \code{\link{flowBasis}} and \code{\link{flowPca}}
}

