### flowFDA: Flow Cytometry Functional Data Analysis 
Authors: Lieven Clement, Olivier Thas, Nico Boon

flowFDA is an [R](http://www.r-project.org) package for analysing flow cytometry
experiments with model based clustering, functional principal component
analysis, functional discriminant analysis and to compare multivariate flow
cytometry fingerprints accross treatments (De Roy, K., Clement, L., Thas, O., Wang, Y., and Boon, N. (2012). Flow cytometry for fast microbial community fingerprinting. Water Research, 46 (3), 907-919.). It will be released on [R/Bioconductor](https://www.bioconductor.org/) in the future. The package is developed at [Ghent University](http://www.ugent.be).


#### Installation

Install flowFDA from its
[GitHub repository](https://github.com/lievenclement/flowFDA). You first need to
install the [devtools](https://cran.r-project.org/package=devtools).

```r
install.packages("devtools")
```

Then install flowFDAE using the `install_github` function in the
[devtools](https://cran.r-project.org/package=devtools) package. (With
`build_vignettes=TRUE`, the vignettes will be built and installed.)

```r
library(devtools)
install_github("lievenclement/flowFDA", build_vignettes=TRUE)
```

### Copyright
Copyright (C) 2016 Lieven Clement, Olivier Thas, Nico Boon.

### Licenses
The flowFDA package as a whole is distributed under
[GPL-3 (GNU General Public License version 3)](http://www.gnu.org/licenses/gpl-3.0.en.html).

