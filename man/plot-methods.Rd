\docType{methods}
\name{plot}
\alias{plot}
\alias{plot,flowBasis-method}
\alias{plot,flowDa-method}
\alias{plot,flowPca-method}
\title{Plot methods
Generate plots for  \code{\link{flowBasis-class}}, \code{\link{flowPca-class}} or  \code{\link{flowDa-class}} objects object}
\arguments{
  \item{x}{\code{\link{flowBasis-class}},
  \code{\link{flowPca-class}} or \code{\link{flowDa-class}}
  objects object.}

  \item{y}{obsolete used for compatibility with S3 plot
  function.}

  \item{samples}{the flows that will be used for
  constructing the two dimensional projections of the
  multivariate distribution. samples=1:5 will plot the
  average projection for the first 5 samples. To be used
  for flowBasis method.}

  \item{pow}{parameter for the coloring cfr. the function
  smoothScatter.}

  \item{ask}{if TRUE the user will be asked to press a
  button before the next plot is displayed.}

  \item{L}{a matrix with a contrasts between flows for
  constructing contrast plots. The different contrasts are
  organized in the columns.}

  \item{contour}{can be used to add contours on the plots.
  The default is FALSE.}

  \item{main=main,...}{additional arguments can be passed
  to customize the plots. e.g. cex, cex.axis, etc.}

  \item{disc}{selects the Principal Components or
  Discriminants that will be visualized}

  \item{nRound}{optional to refine plots (axis labels of
  percentage of variability (discrimination) explained by
  PC (Discriminant)}

  \item{plotType}{different plotTypes can be used.
  plotType="pcaPlot" plots the data in the principal
  component space using standard biplots, plotType="pcaInt"
  produces loading plots for interpretating principal
  components in the original space, plotType="pcaCont"
  produces plots to interpret contrasts in the PCA space,
  plotType="discPlot" produces discriminant biplots,
  plotType="discPlot2" produces stact discriminant plots
  similar to the De Roy et al. (2012) paper,
  plotType="discInt" produces loading plots for
  interpretating discriminants in the original space and
  plotType="discCont" produces plots to interpret contrasts
  in the discriminant space with respect to features of the
  original bivariate projections of the multivariate
  distributions.}

  \item{colRamp}{user defined colorRamp pallette can be
  used when producing interpertation plots for flowPca or
  flowDa objects}

  \item{colorLim}{optional argument for refining
  interpertation plots for flowPca or flowDa objects}

  \item{xAt}{optional to refine plots}

  \item{gLab}{optional to refine plots}

  \item{plotBox}{optional to refine plots}

  \item{groups}{factor with group labels. Optional.}

  \item{cex.Axis}{is used to refine flowDa plots}

  \item{cex.Text}{can be used to refine flowDa plots}
}
\value{
  Noting. Side-effect: plot graphs.
}
\description{
  Plot methods Generate plots for
  \code{\link{flowBasis-class}},
  \code{\link{flowPca-class}} or \code{\link{flowDa-class}}
  objects object
}
\examples{
###########################################
###flowBasis plots
###########################################


if(require(flowFDAExampleData)){
data(fbasis)
#plot of first 6 samples
plot(fbasis,samples=1:6)

#Contrast between control samples and nutrient 3h
L=rep(0,30)
L[1:6]=-1/6
L[25:30]=1/6
plot(fbasis,L=L)
}
#
#
###########################################
###flowDa plots
###########################################

if(require(flowFDAExampleData)){
#Construct fPca object
fPca=flowPca(fbasis)

#Make plots
plot(fPca,disc=1:3)
plot(fPca,plotType="pcaInt",disc=1)
}
#
#
###########################################
###flowDa plots
###########################################
if(require(flowFDAExampleData)){
data(group)
data(fbasis)
#Construct flowDa object
fDa=flowDa(fbasis,groups=group,nPca=.95)
#Make plots, the flowDa class extends flowPca (PCA is performed before discriminant analysis). Hence all plots for flowPca objects can be made for flowDa objects.
plot(fDa,disc=1:3,plotType="pcaPlot")
plot(fDa,disc=1:3)
plot(fDa,disc=1:3,plotType="discPlot2")
plot(fDa,disc=1,plotType="discInt")
}
}

