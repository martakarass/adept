
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


test_that("Test get.x.smoothed, a wrapper for RunningWinSmooth", {

  set.seed(20191021)

  ## Generate data
  x <- sin(seq(0, 2 * pi, length.out = 1000)) + rnorm(1000, sd = 0.2)
  x.smoothed  <- get.x.smoothed(x, 20)
  x.smoothed2 <- RunningWinSmooth(x, W = 20)
  x.smoothed2.na <- which(!is.na(x.smoothed2))

  ## Check that all values from non-NA places do match
  ## between RunningWinSmooth() function and its wrapper, get.x.smoothed()
  expect_true(all(x.smoothed[x.smoothed2.na] == x.smoothed2[x.smoothed2.na]))
  ## Check that there is no NA in the output of get.x.smoothed()
  expect_false(any(is.na(x.smoothed)))
  ## Check that the  output of get.x.smoothed() has its mean unchanged
  expect_equal(mean(x.smoothed), 0.00452722326391469)


  set.seed(20191021)

  ## Generate data
  x <- sin(seq(0, 2 * pi, length.out = 10000)) + rnorm(10000, sd = 0.2)
  x.smoothed  <- get.x.smoothed(x, 99)
  x.smoothed2 <- RunningWinSmooth(x, W = 99)
  x.smoothed2.na <- which(!is.na(x.smoothed2))

  ## Check that all values from non-NA places do match
  ## between RunningWinSmooth() function and its wrapper, get.x.smoothed()
  expect_true(all(x.smoothed[x.smoothed2.na] == x.smoothed2[x.smoothed2.na]))
  ## Check that there is no NA in the output of get.x.smoothed()
  expect_false(any(is.na(x.smoothed)))
  ## Check that the  output of get.x.smoothed() has its mean unchanged
  expect_equal(mean(x.smoothed), 0.001819928)
})


