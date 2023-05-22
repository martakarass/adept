#' Compute sliding covariance like dvmisc::sliding_cov through one convolution.
#'
#' The covariance between two vectors can be written as a single convolution
#' instead of two, as was being done by the dvmisc package.
#'
#' Behavior should be identical to dvmisc's version.
#'
#' @param short Shorter numeric vector to slide over `long`.
#' @param long Numeric vector.
#'
#'
#' @return Numeric vector of sliding covariances.
#'
#' @useDynLib adept
#' @importFrom Rcpp sourceCpp
#'
#' @noRd
#'
slidingCovFast <- function(short, long) {
  n <- length(short)
  len.diff <- length(long) - n
  return(convolveCpp(long, rev(short / (n - 1)))[n:(n + len.diff)])
}


#' Call optimized version of dvmisc::sliding_cor which stores sds.
#'
#' See Rcpp function for more documentation.
#'
#' Assumes mean and sd of `short` are 0 and 1, respectively.
#'
#' @param short Shorter numeric vector to slide over `long`.
#' @param long Numeric vector.
#'
#' @return List of length 2, with $core being a Numeric vector of sliding
#' correlations, the same as what dvmisc::sliding_cor() would output,
#' and $sds being the sliding standard deviations of `long`
#'
#' @useDynLib adept
#' @importFrom Rcpp sourceCpp
#'
#' @noRd
#'
slidingCorStoreSd <- function(short, long) {
  return(slidingCorStoreSdCpp(short, long, 1))
}

#' Call optimized version of dvmisc::sliding_cor.
#'
#' See Rcpp function for more documentation.
#'
#' Assumes mean and sd of `short` are 0 and 1, respectively.
#'
#' @param short Shorter numeric vector to slide over `long`.
#' @param long Numeric vector.
#'
#' @return Numeric vector of sliding correlations, the same as what
#' dvmisc::sliding_cor() would output.
#'
#' @useDynLib adept
#' @importFrom Rcpp sourceCpp
#'
#' @noRd
#'
slidingCor <- function(short, long, sds) {
  return(slidingCorCpp(short, long, 1, sds))
}
