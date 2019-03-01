
[![Travis build status](https://travis-ci.com/martakarass/adept.svg?branch=master)](https://travis-ci.com/martakarass/adept) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/martakarass/adept?branch=master&svg=true)](https://ci.appveyor.com/project/martakarass/adept) [![Coverage status](https://codecov.io/gh/martakarass/adept/branch/master/graph/badge.svg)](https://codecov.io/github/martakarass/adept?branch=master)

<!-- README.md is generated from README.Rmd. Please edit that file -->
adept
=====

`adept` package implements ADaptive Empirical Pattern Transformation (ADEPT) (\[1\]) **method for pattern segmentation from a time-series `x`**. ADEPT was designed for optimal use in performing fast, accurate walking strides segmentation from high-density data collected from a wearable accelerometer worn during continuous walking activity.

Installation
------------

You can install the released version of adept from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("adept")
```

Vignettes
---------

Vignettes are available to better explain package methods functionality.

1.  Vignette [Introduction to adept package](https://martakarass.github.io/adept/articles/adept-intro.html) intends to introduce a reader to the ADEPT method and demonstrate the usage of the `segmentPattern` function which implements ADEPT method. Here, we focus on illustrating `segmentPattern` functionality via simulated data examples.

2.  Vignette [Walking strides segmentation with adept](https://martakarass.github.io/adept/articles/adept-strides-segmentation.html) provides an example of walking stride segmentation from subsecond accelerometry data with `adept` package. We demonstrate that ADEPT method can be used to perform automatic and precise walking stride segmentation from data collected during a combination of running, walking and resting exercise. We demonstrate how to segment data

    -   with the use of stride templates that were pre-computed based on data from an external study,
    -   by deriving new stride templates in a semi-manual manner.
