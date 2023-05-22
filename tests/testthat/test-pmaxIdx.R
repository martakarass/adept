context("Testing pmaxIdx function.")

test_that("pmaxIdx agrees with pmax", {
  input <- list(c(1,2,3,4,5,6),
               c(0,3,0,3,0,3),
               c(0,3,0,6,0,9))
  expected <- c(1,3,3,6,5,9)
  expect_equal(pmaxIdxCpp(input)$pmax, do.call(pmax, input))
  expect_equal(pmaxIdxCpp(input)$pmax, expected)
})

test_that("pmaxIdx agrees with old template index logic", {
  input <- list(c(1,2,3,4,5,6),
               c(0,3,0,3,0,3),
               c(0,3,0,6,0,9))
  expected <- c(1,2,1,3,1,3)
  expect_equal(pmaxIdxCpp(input)$idx, c(max.col(t(do.call(rbind, input)), ties.method = "first")))
  expect_equal(pmaxIdxCpp(input)$idx, expected)
})
