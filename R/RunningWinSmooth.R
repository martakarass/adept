
#' @title
#' Fast Computation of Moving Window Average
#'
#' @description
#' Computes moving window average of a vector via Fast Fourier Transform.
#' The tails of the output vector where the moving window is undefined are filled with \code{NA}.
#'
#' @param x A numeric vector.
#' @param W A width of a moving window given in time (seconds).
#' @param x.fs Frequency of \code{x} expressed in number of observations collected
#' per second. Defaults to \code{1}.
#'
#' @details
#' Frequency of \code{x} and a width of a moving window given in time (seconds) determines
#' \code{W.vl = W * x.fs}, a width of a moving window given in vector length.
#' If \code{W.vl} is  \eqn{W.vl < 3} then an error is thrown. If \code{W.vl} is an even number, then the
#' value \code{W.vl-1} is silently used as a width of a moving window instead.
#'
#' Implementation uses convolution computed via Fast Fourier Transform,
#' which is expected to reduce computational time, especially
#' for a long \code{x} vector, compared to conventional convolution computation.
#'
#' @return A numeric vector of moving window average.
#'
#' @importFrom stats convolve
#'
#' @export
#'
#' @examples
#' ## trivial example f(x) = x
#' N <- 100
#' W  <- 20
#' x <- 1:N
#' x.smoothed <- RunningWinSmooth(x, W)
#' \dontrun{
#' plot(x, type = "l")
#' points(RunningWinSmooth(x, W), col = "red")
#' }
#' ## f(x) = sin(x) + noise
#' N <-  1000
#' W  <- 100
#' x <- sin(seq(0, 4 * pi, length.out = N)) + rnorm(N, sd = 0.1)
#' x.smoothed <- RunningWinSmooth(x, W)
#' \dontrun{
#' plot(x, type = "l")
#' points(RunningWinSmooth(x, W), col = "red")
#' }
#'
RunningWinSmooth <- function(x, W, x.fs = 1){

  W.vl <- round(W * x.fs)
  if (W.vl < 3) stop("W must be not smaller than 3 vector indices. Define wider averaging window length")

  ## Replace W with closest odd integer no larger than W
  W.vl <-  W.vl + (W.vl %% 2) - 1

  ## Comoute moving average via convolution of signal and a fixed value vector
  N <- length(x)
  win <- rep(1/W.vl, W.vl)
  win <- append(win, rep(0, N - W.vl))
  x.out0 <- convolve(x, win)
  x.out0 <- x.out0[1:N]

  ## Redefine head and tail of a signal
  W.wing <- floor(W.vl/2)
  x.out.head <- rep(NA, W.wing)
  x.out.tail <- rep(NA, W.wing)

  x.out <- c(x.out.head,
             x.out0[1:(N - 2 * W.wing)],
             x.out.tail)
  return(x.out)
}






#' Wrapper function for \code{RunningWinSmooth}
#'
#' Wrapper function for \code{RunningWinSmooth}. Replaces \code{NA} values with
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
#' @export
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
      x.smoothed.i <- RunningWinSmooth(x = x.split.i, W = W, x.fs = x.fs)
      x.smoothed.i <- x.smoothed.fill.NA(x.smoothed.i)
      return(x.smoothed.i)
    })
    x.smoothed <- unlist(x.smoothed.list)

  } else {
    x.smoothed <- RunningWinSmooth(x = x, W = W, x.fs = x.fs)
    x.smoothed <- x.smoothed.fill.NA(x.smoothed)
  }

  return(x.smoothed)
}


