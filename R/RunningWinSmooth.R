
#' @title
#' Fast Computation of Moving Window Average
#'
#' @description
#' Compute moving window average of a time-series \code{x}.
#' Head and tail of the output vector where the moving window is undefined are filled with \code{NA}.
#'
#' @param x A numeric vector. A time-series for which moving window average is computed.
#' @param W A numeric scalar. A length of a moving window, expressed in time (seconds).
#' @param x.fs Frequency of time-series \code{x}, expressed in a number of observations
#' per second. Defaults to \code{1}.
#'
#' @details
#' Time-series frequency \code{x.fs} and a length  of a moving window (expressed in time) \code{W}
#' together determine
#' \code{W.vl = round(W * x.fs)}, a length of a moving window expressed in \code{x} vector length.
#' Note: \code{W.vl} must be equal or greater than \code{3}.
#' \itemize{
#'   \item If \code{W.vl < 3} then an error is thrown.
#'   \item If \code{W.vl} is an even number then \code{(W.vl-1)} value is silently
#'   used instead as a length of a moving window expressed in \code{x} vector length.
#' }
#'
#' @return A numeric vector of a moving window average.
#'
#' @importFrom stats convolve
#'
#' @export
#'
#' @examples
#' ## Time-series defined as a function f(x) = x
#' N <- 100
#' W  <- 20
#' x <- 1:N
#' x.smoothed <- runningWinSmooth(x, W)
#' \dontrun{
#' plot(x, type = "l")
#' points(x.smoothed, col = "red")
#' }
#'
#' ## Time-series defined as a function f(x) = sin(x) + noise
#' N <-  1000
#' W  <- 100
#' x <- sin(seq(0, 4 * pi, length.out = N)) + rnorm(N, sd = 0.1)
#' x.smoothed <- runningWinSmooth(x, W)
#' \dontrun{
#' plot(x, type = "l")
#' points(x.smoothed, col = "red")
#' }
#'
runningWinSmooth <- function(x, W, x.fs = 1){

  ## Check function arguments for correctness
  if ((!(is.vector(x))) || (!(is.vector(x)))) stop("x must be a numeric vector.")
  if (!is.numeric(W)) stop("W must be a numeric scalar")
  if (length(W) > 1) stop("W must be a numeric scalar (1-element numeric vector)")
  if (!is.numeric(x.fs)) stop("x.fs must be a numeric scalar")
  if (length(x.fs) > 1) stop("x.fs must be a numeric scalar (1-element numeric vector)")

  W.vl <- round(W * x.fs)
  if (W.vl < 3) stop("W.vl (refer to function's details description) must not be smaller than 3 vector indices. Define wider W averaging window length")
  if (W.vl >= length(x)) stop("W.vl (refer to function's details description) must be less than x vector length. Define narrower W averaging window length")

  ## Replace W with closest odd integer no larger than W
  W.vl <-  W.vl + (W.vl %% 2) - 1

  ## Compute moving average via convolution of signal and a fixed value vector
  N <- length(x)
  win <- rep(1/W.vl, W.vl)
  win <- append(win, rep(0, N - W.vl))
  x.out0 <- convolve(x, win)
  x.out0 <- x.out0[1:N]

  ## Replace head and tail of a signal with NA's
  W.wing <- floor(W.vl/2)
  x.out.head <- rep(NA, W.wing)
  x.out.tail <- rep(NA, W.wing)

  x.out <- c(x.out.head,
             x.out0[1:(N - 2 * W.wing)],
             x.out.tail)
  return(x.out)
}




#' Wrapper function for \code{runningWinSmooth}
#'
#' Wrapper function for \code{runningWinSmooth}. Replaces \code{NA} values with
#' which appear in head and tail of smoothed signal as a result of MA not defined
#' for first and last \code{k} elements of a signal. The \code{NA} values are replaced
#' with sample means of subsequent/procceeding non-\code{NA} values that appear
#' in smoothed signal.
#'
#'
#' @param x Numeric vector to be smoothed.
#' @param W A width of a moving window given in time (seconds).
#' @param x.fs Frequency of \code{x} expressed in number of observations collected
#' per second. Defaults to \code{1}.
#' @param NA.repl.source.k Scalar for number of subsequent/procceeding non-\code{NA} values that appear
#' in smoothed signal and are used to provide sample mean replacement for \code{NA} values.
#'
#' @return Smoothed signal with  \code{NA} values replaced with sample means of
#' subsequent/procceeding non-\code{NA} values that appear
#' in smoothed signal.
#'
#' @examples
#' set.seed(20191021)
#' x <- sin(seq(0, 2 * pi, length.out = 1000)) + rnorm(1000, sd = 0.2)
#' x.smoothed <- get.x.smoothed(x, 20)
#' plot(x, type = "l", col = "grey")
#' lines(x.smoothed, col = "blue")
#'
#' @noRd
#'
get.x.smoothed <- function(x, W, x.fs = 1, NA.repl.source.k = 4,
                           x.cut.vl = 100 * 60 * 60 * 24){

  W.vl <- W * x.fs

  ## Nested fuction responsible for filling the NAs
  x.smoothed.fill.NA <- function(x.smoothed){
    ## Replace NA's in head/tail of smoothed signal with some neutral average flat line
    ## Vector length of replacement NA's area
    NA.vl <- floor((W.vl + (W.vl %% 2) - 1)/2)
    ## Replace NAs in vector's head
    repl.head.idx        <- 1:NA.vl
    repl.head.source.idx <- 1:(NA.repl.source.k * NA.vl)
    x.smoothed[repl.head.idx] <- mean(x.smoothed[repl.head.source.idx], na.rm = TRUE)
    ## Replace NAs in vector's tail
    repl.tail.idx        <- length(x.smoothed) - ((NA.vl - 1):0)
    repl.tail.source.idx <- length(x.smoothed) - ((NA.repl.source.k * NA.vl - 1):0)
    x.smoothed[repl.tail.idx] <- mean(x.smoothed[repl.tail.source.idx], na.rm = TRUE)
    return(x.smoothed)
  }

  ## Decide whether to smooth on x parts or not
  x.vl <- length(x)
  if (x.vl > x.cut.vl){
    rn <- ceiling(x.vl / x.cut.vl)
    x.split.idx <- unlist(lapply(1:rn, function(i) rep(i, x.cut.vl)))[1:x.vl]
    if ((sum(x.split.idx == rn) < 100 * 60) & (rn > 1)){
      x.split.idx[x.split.idx == rn] <- rn - 1
    }
    x.split <- split(x, x.split.idx)
    ## Apply to each part of x split
    x.smoothed.list <- lapply(x.split, function(x.split.i){
      x.smoothed.i <- runningWinSmooth(x = x.split.i, W = W, x.fs = x.fs)
      x.smoothed.i <- x.smoothed.fill.NA(x.smoothed.i)
      return(x.smoothed.i)
    })
    x.smoothed <- unlist(x.smoothed.list)

  } else {
    x.smoothed <- runningWinSmooth(x = x, W = W, x.fs = x.fs)
    x.smoothed <- x.smoothed.fill.NA(x.smoothed)
  }

  return(x.smoothed)
}

# get.x.smoothed<- function(x, W, x.fs = 1, NA.repl.surce.k = 4){
#
#   W.vl <- W * x.fs
#   x.smoothed <- runningWinSmooth(x = x, W = W, x.fs = x.fs)
#   ## Replace NA's in head/tail of smoothed signal with some neutral average flat line
#   ## Vector length of replacement NA's area
#
#   NA.vl <- floor((W.vl + (W.vl %% 2) - 1)/2)
#   ## Replace NAs in vector's head
#   repl.head.idx        <- 1:NA.vl
#   repl.head.source.idx <- 1:(NA.repl.surce.k * NA.vl)
#   x.smoothed[repl.head.idx] <- mean(x.smoothed[repl.head.source.idx], na.rm = TRUE)
#
#   ## Replace NAs in vector's tail
#   repl.tail.idx        <- length(x.smoothed) - ((NA.vl - 1):0)
#   repl.tail.source.idx <- length(x.smoothed) - ((NA.repl.surce.k * NA.vl - 1):0)
#   x.smoothed[repl.tail.idx] <- mean(x.smoothed[repl.tail.source.idx], na.rm = TRUE)
#
#   return(x.smoothed)
#
# }


