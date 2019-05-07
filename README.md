
[![Travis build
status](https://travis-ci.com/martakarass/adept.svg?branch=master)](https://travis-ci.com/martakarass/adept)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/martakarass/adept?branch=master&svg=true)](https://ci.appveyor.com/project/martakarass/adept)
[![Coverage
status](https://codecov.io/gh/martakarass/adept/branch/master/graph/badge.svg)](https://codecov.io/github/martakarass/adept?branch=master)

<!-- README.md is generated from README.Rmd. Please edit that file -->

The `adept` package implements ADaptive Empirical Pattern Transformation
(ADEPT) method\[1\] for pattern segmentation from a time-series. ADEPT
is designed to perform fast, accurate walking strides segmentation from
high-density data collected with a wearable accelerometer during
walking.

### Installation

Install `adept` package from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("martakarass/adept")
```

### Example 1

Simulate a time-series `x`; assume that `x` is collected frequency of
100 Hz, there is one pattern of a fixed duration 1.01 seconds present in
`x`, and there is no noise in collected data.

``` r
library(adept)

x <- cos(seq(0, 2 * pi * 10, length.out = 1001))
true.pattern <- x[1:101]

par(mfrow = c(1,2), cex = 1)
plot(true.pattern, type = "l", xlab = "", ylab = "", main = "Pattern")
plot(x, type = "l", xlab = "", ylab = "", main = "Time-series x")
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

Segment pattern from data. Assume perfect template is available. Use
grid of potential pattern durations of {0.9, 1, 1.1} seconds.

``` r
segmentPattern(
  x = x,
  x.fs = 100,
  template = true.pattern,
  pattern.dur.seq = c(0.9, 1, 1.1),
  similarity.measure = "cor",
  compute.template.idx = TRUE)
#>    tau_i T_i     sim_i template_i
#> 1      1 100 0.9994793          1
#> 2    101 100 0.9994793          1
#> 3    201 100 0.9994793          1
#> 4    302 100 0.9994793          1
#> 5    401 100 0.9994793          1
#> 6    501 100 0.9994793          1
#> 7    601 100 0.9994793          1
#> 8    701 100 0.9994793          1
#> 9    801 100 0.9994793          1
#> 10   901 100 0.9994793          1
```

Assume grid of potential pattern durations of {0.9, 1.01, 1.1} seconds
(the grid is now “perfect” in a sense it contains duration of the true
pattern used in data simulation).

``` r
segmentPattern(
  x = x,
  x.fs = 100,
  template = true.pattern,
  pattern.dur.seq = c(0.9, 1.01, 1.1),
  similarity.measure = "cor",
  compute.template.idx = TRUE)
#>    tau_i T_i sim_i template_i
#> 1      1 101     1          1
#> 2    101 101     1          1
#> 3    201 101     1          1
#> 4    301 101     1          1
#> 5    401 101     1          1
#> 6    501 101     1          1
#> 7    601 101     1          1
#> 8    701 101     1          1
#> 9    801 101     1          1
#> 10   901 101     1          1
```

### Example 2

Simulate a time-series `x`; assume that time-series is collected
frequency of 100 Hz, there are two patterns of various duration present
in `x`, and there is noise in collected data.

Then, generate `x2` as a noisy version of `x`.

``` r
true.pattern.1 <- cos(seq(0, 2 * pi, length.out = 200))
true.pattern.2 <- true.pattern.1
true.pattern.2[70:130] <- 2 * true.pattern.2[min(70:130)] + abs(true.pattern.2[70:130])
x <- numeric()
for (vl in seq(70, 130, by = 10)){
  true.pattern.1.s <- approx(
    seq(0, 1, length.out = 200), 
    true.pattern.1, xout = seq(0, 1, length.out = vl))$y
  true.pattern.2.s <- approx(
    seq(0, 1, length.out = 200), 
    true.pattern.2, xout = seq(0, 1, length.out = vl))$y
  x <- c(x, true.pattern.1.s[-1], true.pattern.2.s[-1])
  if (vl == 70) x <- c(true.pattern.1.s[1], x)
}
set.seed(1)
x2 <- x + rnorm(length(x), sd = 0.5)

par(mfrow = c(1,2), cex = 1)
plot(true.pattern.1, type = "l", xlab = "", ylab = "", main = "Pattern 1")
plot(true.pattern.2, type = "l", xlab = "", ylab = "", main = "Pattern 2")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r
par(mfrow = c(1,1), cex = 1)
plot(x, type = "l", xlab = "", ylab = "", main = "Time-series x")
```

<img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />

``` r
plot(x2, type = "l", xlab = "", ylab = "", main = "Time-series x2")
```

<img src="man/figures/README-unnamed-chunk-5-3.png" width="100%" />

Assume perfect grid of potential pattern duration, {0.7, 0.8, 0.9, 1.0,
1.1, 1.2, 1.3} seconds. Segment `x`.

``` r
segmentPattern(
  x = x,
  x.fs = 100,
  template = list(true.pattern.1, true.pattern.2),
  pattern.dur.seq = seq(70, 130, by = 10) * 0.01,
  similarity.measure = "cor",
  compute.template.idx = TRUE)
#>    tau_i T_i sim_i template_i
#> 1      1  70     1          1
#> 2     70  70     1          2
#> 3    139  80     1          1
#> 4    218  80     1          2
#> 5    297  90     1          1
#> 6    386  90     1          2
#> 7    475 100     1          1
#> 8    574 100     1          2
#> 9    673 110     1          1
#> 10   782 110     1          2
#> 11   891 120     1          1
#> 12  1010 120     1          2
#> 13  1129 130     1          1
#> 14  1258 130     1          2
```

Segment `x2`.

``` r
segmentPattern(
  x = x2,
  x.fs = 100,
  template = list(true.pattern.1, true.pattern.2),
  pattern.dur.seq = seq(70, 130, by = 10) * 0.01,
  similarity.measure = "cor",
  compute.template.idx = TRUE)
#>    tau_i T_i     sim_i template_i
#> 1      1  70 0.8585451          1
#> 2    138  80 0.7624002          1
#> 3    218  80 0.7025577          2
#> 4    297  90 0.8500864          1
#> 5    390  80 0.6931671          2
#> 6    469 110 0.8286013          1
#> 7    579  90 0.6373846          2
#> 8    668 120 0.8027177          1
#> 9    787 100 0.6666713          2
#> 10   888 130 0.7894766          1
#> 11  1017 110 0.6599280          1
#> 12  1129 130 0.7938183          1
#> 13  1267 120 0.7655408          2
```

Segment `x2`; smooth data time-series before similarity matrix
computation. Assume more dense grid of potential pattern duration.

``` r
par(mfrow = c(1,1), cex = 1)
plot(windowSmooth(x = x2, x.fs = 100, W = 0.1), 
     type = "l", xlab = "", ylab = "", main = "Time-series x2 smoothed")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

``` r
segmentPattern(
  x = x2,
  x.fs = 100,
  template = list(true.pattern.1, true.pattern.2),
  pattern.dur.seq = 70:130 * 0.01,
  similarity.measure = "cor",
  x.adept.ma.W = 0.1,
  compute.template.idx = TRUE)
#>    tau_i T_i     sim_i template_i
#> 1      1  70 0.9865778          1
#> 2     70  70 0.9533684          2
#> 3    139  79 0.9683054          1
#> 4    217  80 0.9748040          2
#> 5    296  94 0.9802473          1
#> 6    391  82 0.9462213          2
#> 7    472 106 0.9855837          1
#> 8    578  93 0.9608881          2
#> 9    670 115 0.9887225          1
#> 10   784 107 0.9562694          2
#> 11   896 113 0.9734575          1
#> 12  1008 127 0.9703118          1
#> 13  1134 116 0.9606235          1
#> 14  1266 122 0.9593345          2
```

### Vignettes

Vignettes are available to better demonstrate package methods usgae.

1.  Vignette [Introduction to adept
    package](https://martakarass.github.io/adept/articles/adept-intro.html)
    introduces ADEPT algorithm and demonstrates the usage of
    `segmentPattern` function which implements ADEPT approach. Here, we
    focus on examples with simulated data.

2.  Vignette [Walking strides segmentation with
    adept](https://martakarass.github.io/adept/articles/adept-strides-segmentation.html)
    provides an example of segmentation of walking strides (two
    consecutive steps) in sub-second accelerometry data with `adept`
    package. The exemplary dataset is a part of the `adeptdata` package.
    We demonstrate that ADEPT can be used to perform automatic and
    precise walking stride segmentation from data collected during a
    combination of running, walking and resting exercises. We introduce
    how to segment data:
    
    1.  with the use of stride templates that were pre-computed based on
        data from an external study (attached to `adeptdata` package),
    2.  by deriving new stride templates in a semi-manual manner.

References

1.  Karas, M., Straczkiewicz, M., Fadel, W., Harezlak, J., Crainiceanu,
    C., Urbanek, J.K. *Adaptive empirical pattern transformation (ADEPT)
    with application to walking stride segmentation*, Submitted to
    *Biostatistics*, 2018.
