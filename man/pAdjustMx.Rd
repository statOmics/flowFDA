\name{pAdjustMx}
\alias{pAdjustMx}
\title{Function that extends p.adjust to matrices of p-values}
\usage{
  pAdjustMx(pvals, method = "holm")
}
\arguments{
  \item{pvals}{matrix with p-values}

  \item{method}{method to correct for multiple testing,
  c("holm", "hochberg", "hommel", "bonferroni", "BH",
  "BY","fdr", "none"). Details can be found in
  stats::p.adjust.}
}
\value{
  matrix with adjusted pvalues
}
\description{
  Function that extends p.adjust to matrices of p-values
}
\examples{
if(require(flowFDAExampleData)){
data(fDa)
pAdjustMx(getMpc(fDa)$p)
}
}

