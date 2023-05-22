context("Testing sliding cov and cor functions.")

test_that("slidingCovFast agrees with dvmisc::sliding_cov on short vectors.", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(1.5, 1, 2.1)
  expect_equal(dvmisc::sliding_cov(short, long),
               slidingCovFast(short, long))
})

test_that("slidingCovFast agrees with dvmisc::sliding_cov on long vectors.", {
  set.seed(1)
  long <- runif(1000)
  short <- runif(500)
  expect_equal(dvmisc::sliding_cov(short, long),
               slidingCovFast(short, long))
})

test_that("slidingCorStoreSd agrees with dvmisc::sliding_cor.", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(1.5, 1, 2.1)
  expect_equal(dvmisc::sliding_cor(short, long),
               slidingCorStoreSd(short, long)$cor)
})

test_that("slidingCorStoreSd returns NA for constant sections.", {
  long <- c(1, 5, 2, 7, 7, 7)
  short <- c(1.5, 1, 2.1)
  expect_equal(c(dvmisc::sliding_cor(short, long)[1:3], NA),
               slidingCorStoreSd(short, long)$cor)
})

test_that("slidingCor agrees with slidingCorStoreSd.", {
  long <- c(1, 5, 2, 7, 7, 7)
  short <- c(1.5, 1, 2.1)
  expected <- slidingCorStoreSd(short, long)
  expect_equal(expected$cor,
               slidingCor(short, long, expected$sds))
})
