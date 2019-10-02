[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![license](https://img.shields.io/badge/license-GPL3-lightgrey.svg)](https://choosealicense.com/)
[![DOI](https://zenodo.org/badge/173801188.svg)](https://zenodo.org/badge/latestdoi/173801188)

rcr
===

rcr is an R package to calculate a basic backwater profile for open channels using the standard step method.

The package is currently configured for backwater calculations with the subcritical assumption in open channels only, supercritical and mixed regimes, and hydrualic structures, will be added in future versions.

Helpful functions in the package include:

-   *compute\_flow\_profiles* computes flow profiles for a given geometry for a vector of flow values
-   *bankfull\_estimator* applies Manning's equation to estimate flow for a given cross-section, slope and water surface level
-   *flow\_area* calculates the flow area for a set of stations in a cross-section
-   *normal\_depth* calculates the normal depth of a cross-section

Installation
------------

You can install rcr from github with:

``` r
# install.packages("devtools")
devtools::install_github("rchlumsk/rcr")
```

rcr Wishlist
------------

Any issues or feature requests can be submitted on Github via the Issues tab. You may also submit feature requests directly to Robert Chlumsky (<rchlumsk@uwaterloo.ca>) via email.

Examples
--------

See the vignettes in the package for examples of using the package functions using the sample data provided with the package.
