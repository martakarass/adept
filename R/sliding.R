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
  return(convolveCpp(long, rev(short / (n - 1) - sum(short) / n / (n - 1)))[n:(n + len.diff)])
}

#' Call optimized version of dvmisc::sliding_cor.
#'
#' The original version had one completely redundant sum, and was recomputing
#' parts of sums redundantly as well. This calls a custom Rcpp version that
#' does neither.
#'
#' Behavior should be identical to dvmisc's version.
#'
#' @param short Shorter numeric vector to slide over `long`.
#' @param long Numeric vector.
#'
#'
#' @return Numeric vector of sliding correlations.
#'
#' @useDynLib adept
#' @importFrom Rcpp sourceCpp
#'
#' @noRd
#'
slidingCorFast <- function(short, long) {
  return(slidingCorCpp(short, long, sd(short)))
}
