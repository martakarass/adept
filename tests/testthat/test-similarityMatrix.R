

context("Testing similarityMatrix function.")

test_that("Testing nothing has change in similarityMatrix() output for a case with 1 template, 3 scales.", {

  set.seed(20181018)

  ## Simulate data
  x0 <- sin(seq(0, 2 * pi * 10, length.out = 10000))
  x  <- x0 + rnorm(10000, sd = 0.1)
  template <- list(x0[1:500])
  template.vl <- c(400, 500, 600)
  template.scaled <- scaleTemplate(template, template.vl)


  ## Test for similarity.measure = "cov"
  out <- similarityMatrix(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               -0.000104255166378453)

  ## Test for similarity.measure = "cor"
  out <- similarityMatrix(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
               -0.000610391199176645)

})





test_that("Testing nothing has change in similarityMatrix() output for a case with 3 templates, 3 scales.", {

  set.seed(20181018)

  ## Simulate data
  x0 <- sin(seq(0, 2 * pi * 10, length.out = 10000))
  x  <- x0 + rnorm(10000, sd = 0.1)
  template <- list(x0[1:500], x0[51:550], x0[51:400])
  template.vl <- c(400, 500, 600)
  template.scaled <- scaleTemplate(template, template.vl)

  ## Test for similarity.measure = "cov"
  out <- similarityMatrix(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               0.240442503501679)

  ## Test for similarity.measure = "cor"
  out <- similarityMatrix(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
               0.422636316622436)
})




test_that("Testing nothing has change in similarityMatrix() output for a case with 3 templates, 3 scales.", {

  set.seed(20181018)

  ## Simulate data
  x0 <- sin(seq(0, 2 * pi * 10, length.out = 10000))
  x  <- x0 + rnorm(10000, sd = 0.1)
  template <- list(x0[1:500], x0[51:550], x0[51:400])
  template.vl <- c(400, 500, 600)
  template.scaled <- scaleTemplate(template, template.vl)

  ## Test for similarity.measure = "cov"
  out <- similarityMatrix(x, template.scaled, "cov")
  expect_equal(mean(out, na.rm = TRUE),
               0.240442503501679)

  ## Test for similarity.measure = "cor"
  out <- similarityMatrix(x, template.scaled, "cor")
  expect_equal(mean(out, na.rm = TRUE),
               0.422636316622436)
})



test_that("Testing nothing has change in templateIdxMatrix() output", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  # s.grid <- sample(60:120, size = 5)
  s.grid <- c(76L, 82L, 93L, 112L, 71L)
  template1 <- cos(seq(0, 2 * pi, length.out = 200))
  template2 <- c(rev(seq(-1, 1, length.out = 100)), seq(-1, 1, length.out = 100))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    ## Add piece from template1
    templ0 <- approx(seq(0, 1, length.out = 200), template1, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0) x <- x[-length(x)]
    x <- c(x, templ0)
    ## Add piece from template2
    templ0 <- approx(seq(0, 1, length.out = 200), template2, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0) x <- x[-length(x)]
    x <- c(x, templ0)
  }

  ## Run ADEPT
  template <- list(template1, template2)
  template.vl <- sort(s.grid)
  template.scaled <- scaleTemplate(template, template.vl)
  similarity.measure <- "cor"
  out <- templateIdxMatrix(x, template.scaled, similarity.measure)

  ## Test table of counts is the same
  res <- table(out)
  res.exp <- structure(c(`1` = 1849L, `2` = 2017L), .Dim = 2L, .Dimnames = list(
    out = c("1", "2")), class = "table")
  expect_equal(res, res.exp)

  ## Test diff same
  res <- apply(out, 1, function(row) mean(diff(row), na.rm = TRUE))
  res.exp <- c(0.00126903553299492, 0.00127713920817369, 0, 0, 0.00133868808567604)
  expect_equal(res, res.exp)
})


