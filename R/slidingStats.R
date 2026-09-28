# Fast native moving-window statistics used internally by adept.

sliding_cov <- function(short, long) {
  if (!is.numeric(short) || !is.atomic(short) || !is.null(dim(short)) ||
      !is.numeric(long) || !is.atomic(long) || !is.null(dim(long))) {
    stop("short and long must be numeric vectors")
  }
  short <- as.double(short)
  long <- as.double(long)
  if (length(short) < 2L) stop("short must contain at least two observations")
  if (length(long) < length(short)) {
    stop("long must be at least as long as short")
  }
  .Call("adept_sliding_cov", short, long, PACKAGE = "adept")
}

sliding_cor <- function(short, long) {
  if (!is.numeric(short) || !is.atomic(short) || !is.null(dim(short)) ||
      !is.numeric(long) || !is.atomic(long) || !is.null(dim(long))) {
    stop("short and long must be numeric vectors")
  }
  short <- as.double(short)
  long <- as.double(long)
  if (length(short) < 2L) stop("short must contain at least two observations")
  if (length(long) < length(short)) {
    stop("long must be at least as long as short")
  }
  .Call("adept_sliding_cor", short, long, PACKAGE = "adept")
}

moving_mean <- function(x, window) {
  if (!is.numeric(x) || !is.atomic(x) || !is.null(dim(x))) {
    stop("x must be a numeric vector")
  }
  if (!is.numeric(window) || length(window) != 1L || is.na(window) ||
      !is.finite(window) || window < 1) {
    stop("window must be a positive integer")
  }
  if (window > length(x)) stop("window must not exceed the length of x")
  if (window != floor(window)) stop("window must be a positive integer")
  x <- as.double(x)
  window <- as.integer(window)
  .Call("adept_moving_mean", x, window, PACKAGE = "adept")
}
