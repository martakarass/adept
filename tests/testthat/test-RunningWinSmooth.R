
context("Testing RunningWinSmooth function.")

test_that("Smoothing of f(x) = x", {

  N <- 100
  W  <- 20
  x <- 1:N
  x.smoothed <- RunningWinSmooth(x, W)

  expect_equal(sum(is.na(x.smoothed)), 18)
  expect_equal(max(cumsum(x.smoothed[!is.na(x.smoothed)])), 4141)
})


test_that("Smoothing of f(x) = sin(x) + noise", {

  set.seed(20180927)

  N <- 1000
  W  <- 200
  x <- sin(seq(0, 4 * pi, length.out = N)) + rnorm(N, sd = 0.1)
  x.smoothed <- RunningWinSmooth(x, W)

  expect_equal(mean(x.smoothed, na.rm = TRUE), 0.0037456166280208)
})


test_that("Smoothing of f(x) = sin(x) + noise, for long signal f(x)", {

  set.seed(20180927)

  N <- 100000
  W  <- 100
  x <- sin(seq(0, 4 * pi, length.out = N)) + rnorm(N, sd = 0.1)
  x.smoothed <- RunningWinSmooth(x, W)

  expect_equal(mean(x.smoothed, na.rm = TRUE), 0.000299133646333034)
})
