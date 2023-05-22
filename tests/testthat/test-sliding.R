context("Testing sliding cov and cor functions.")

test_that("slidingCovFast agrees with dvmisc::sliding_cov", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(1.5, 1, 2.1)
  expect_equal(dvmisc::sliding_cov(short, long),
               slidingCovFast(short, long))
})

test_that("slidingCorFast agrees with dvmisc::sliding_cor", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(1.5, 1, 2.1)
  expect_equal(dvmisc::sliding_cor(short, long),
               slidingCorFast(short, long))
})
