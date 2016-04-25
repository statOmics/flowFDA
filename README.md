### flowFDA: Flow Cytometry Functional Data Analysis 
Authors: Lieven Clement, Olivier Thas, Nico Boon

flowFDA is an [R](http://www.r-project.org) package for analysing flow cytometry
experiments with model based clustering, functional principal component
analysis, functional discriminant analysis and to compare multivariate flow
cytometry fingerprints accross treatments (De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919.). It will be released on [R/Bioconductor](https://www.bioconductor.org/) in the future. The package is developed at [Ghent University](http://www.ugent.be).


#### Installation

Install flowFDA from its
[GitHub repository](https://github.com/lievenclement/flowFDA). You first need to
install the R/Bioconductor packages [Biobase](http://bioconductor.org/packages/release/bioc/html/Biobase.html), [BiocGenerics](http://bioconductor.org/packages/release/bioc/html/BiocGenerics.html) (>= 0.1.6), [flowCore] (http://bioconductor.org/packages/release/bioc/html/flowCore.html), [flowViz](http://bioconductor.org/packages/release/bioc/html/BiocGenerics.html), [flowFP](http://bioconductor.org/packages/release/bioc/html/flowFP.html),
[graphics](https://cran.r-project.org/package=graphics), [grDevices](https://cran.r-project.org/package=grDevices), [methods](https://cran.r-project.org/package=methods), [stats](https://cran.r-project.org/package=stats), [stats4](https://cran.r-project.org/package=stats4), [MASS](https://cran.r-project.org/package=MASS), [multcomp](https://cran.r-project.org/package=multcomp),[mclust](https://cran.r-project.org/package=mclust) and [devtools](https://cran.r-project.org/package=devtools).

```r
source("https://bioconductor.org/biocLite.R")
biocLite(c("flowCore", "flowViz", "flowFP", "MASS", "multcomp","mclust","devtools")
```

Then install flowFDAE using the `install_github` function in the
[devtools](https://cran.r-project.org/package=devtools) package. (With
`build_vignettes=TRUE`, the vignettes will be built and installed.) 
You first need to install the flowFDADataExample package for this purpose

```r
library(devtools)
install_github("lievenclement/flowFDAExampleData")
install_github("lievenclement/flowFDA", build_vignettes=TRUE)
```

### Copyright
Copyright (C) 2016 Lieven Clement, Olivier Thas, Nico Boon.

### Licenses
The flowFDA package as a whole is distributed under
[GPL-3 (GNU General Public License version 3)](http://www.gnu.org/licenses/gpl-3.0.en.html).

