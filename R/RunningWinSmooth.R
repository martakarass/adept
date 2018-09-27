
#' @title
#' Fast Computation of Moving Window Average
#'
#' @description
#' Computes moving window average of a vector via Fast Fourier Transform.
#' The tails of the output vector where the moving window is undefined are filled with \code{NA}.
#'
#' @param x A numeric vector.
#' @param W A width of a moving window given in number of indices. Must be \eqn{\geq 3} and
#' an odd number (see Details).
#'
#' @details
#' If \code{W} is  \eqn{< 3} then an error is thrown. If \code{W} is an even number, then the
#' value \code{W-1} is silently used as a width of a moving window instead.
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
RunningWinSmooth <- function(x, W){

  if (W < 3) stop("W must be not smaller than 3 vector indices. Define wider averaging window length")

  ## Replace W with closest odd integer no larger than W
  W <-  W + (W %% 2) - 1

  ## Comoute moving average via convolution of signal and a fixed value vector
  N <- length(x)
  win <- rep(1/W, W)
  win <- append(win, rep(0, N - W))
  x.out0 <- convolve(x, win)
  x.out0 <- x.out0[1:N]

  ## Redefine head and tail of a signal
  W.wing <- floor(W/2)
  x.out.head <- rep(NA, W.wing)
  x.out.tail <- rep(NA, W.wing)

  x.out <- c(x.out.head,
             x.out0[1:(N - 2 * W.wing)],
             x.out.tail)
  return(x.out)
}
