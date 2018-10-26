
## Note:
## We cannot test on more than 2 cores
## https://stackoverflow.com/questions/41307178/error-processing-vignette-failed-with-diagnostics-4-simultaneous-processes-spa

context("Testing adeptSimilarity function.")

test_that("Testing nothing has change for the result with 1 template only", {

  set.seed(20181018)

  x0 <- sin(seq(0, 2 * pi * 10, length.out = 10000))
  x  <- x0 + rnorm(10000, sd = 0.1)

  template <- list(x0[1:500])
  template.vl <- c(400, 500, 600)
  template.scaled <- scaleTemplate(template, template.vl)

  ## ---------------------------------------------------------------------------
  ## cov

  ## Non-parallel
  out <- adeptSimilarity(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               -0.000104255166378453)

  ## Parallel
  # out2 <- adeptSimilarity(x, template.scaled, "cov",
  #                         run.parallel = TRUE)
  # expect_equal(mean(out2, na.rm = TRUE),
  #              -0.000104255166378453)

  out3 <- adeptSimilarity(x, template.scaled, "cov",
                          run.parallel = TRUE,
                          run.parallel.ncores = 2)
  expect_equal(mean(out3, na.rm = TRUE),
               -0.000104255166378453)


  ## ---------------------------------------------------------------------------
  ## cor

  ## Non-parallel
  out <- adeptSimilarity(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
               -0.000610391199176645)

  ## Parallel
  # out2 <- adeptSimilarity(x, template.scaled, "cor",
  #                         run.parallel = TRUE)
  # expect_equal(mean(out2, na.rm = TRUE),
  #              -0.000610391199176645)

  out3 <- adeptSimilarity(x, template.scaled, "cor",
                          run.parallel = TRUE,
                          run.parallel.ncores = 2)
  expect_equal(mean(out3, na.rm = TRUE),
               -0.000610391199176645)

})





test_that("Testing nothing has change for the result with 3 templates", {

  set.seed(20181018)

  x0 <- sin(seq(0, 2 * pi * 10, length.out = 10000))
  x  <- x0 + rnorm(10000, sd = 0.1)

  template <- list(x0[1:500], x0[51:550], x0[51:400])
  template.vl <- c(400, 500, 600)
  template.scaled <- scaleTemplate(template, template.vl)

  ## ---------------------------------------------------------------------------
  ## cov

  ## Non-parallel
  out <- adeptSimilarity(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               0.240442503501679)

  ## Parallel
  # out2 <- adeptSimilarity(x, template.scaled, "cov",
  #                         run.parallel = TRUE)
  # expect_equal(mean(out2, na.rm = TRUE),
  #              0.240442503501679)

  out3 <- adeptSimilarity(x, template.scaled, "cov",
                          run.parallel = TRUE,
                          run.parallel.ncores = 2)
  expect_equal(mean(out3, na.rm = TRUE),
               0.240442503501679)


  ## ---------------------------------------------------------------------------
  ## cor

  ## Non-parallel
  out <- adeptSimilarity(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
               0.422636316622436)

  ## Parallel
  # out2 <- adeptSimilarity(x, template.scaled, "cor",
  #                         run.parallel = TRUE)
  # expect_equal(mean(out2, na.rm = TRUE),
  #              0.422636316622436)

  out3 <- adeptSimilarity(x, template.scaled, "cor",
                          run.parallel = TRUE,
                          run.parallel.ncores = 2)
  expect_equal(mean(out3, na.rm = TRUE),
               0.422636316622436)

})
