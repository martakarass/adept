

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

  out <- adeptSimilarity(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               -0.000104255166378453)

  ## ---------------------------------------------------------------------------
  ## cor

  out <- adeptSimilarity(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
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

  out <- adeptSimilarity(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               0.240442503501679)


  ## ---------------------------------------------------------------------------
  ## cor
  out <- adeptSimilarity(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
               0.422636316622436)
})
