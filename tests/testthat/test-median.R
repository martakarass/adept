context("Testing Rcpp median function agrees with base R.")


test_that("medianCpp agrees with median without NA.", {
  x = 1:10
  expect_equal(median(x),
               medianCpp(x))
})

test_that("medianCpp agrees with median with NAs.", {
  x = c(NA, 1:10)
  expect_equal(median(x),
               medianCpp(x))
})

test_that("medianCpp agrees with median with NaNs.", {
  x = c(NaN, 1:10)
  expect_equal(median(x),
               medianCpp(x))
})
