
context("Testing scaleTemplate function.")

test_that("Testing correctness of number of elements in output list", {

  set.seed(20180927)

  template    <- list(rnorm(100), rnorm(110), rnorm(120), rnorm(130), rnorm(140))
  template.vl <- seq(100, 200, length.out = 10)
  template.vl <- sort(unique(round(template.vl)))

  res <- scaleTemplate(template, template.vl)

  expect_equal(length(res), 10)
  expect_equal(unlist(lapply(res, length)),
               c(5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L))
})


test_that("Testing correctness of length of elements in output list", {

  set.seed(20180927)

  template    <- list(rnorm(100), rnorm(110), rnorm(120), rnorm(130), rnorm(140))
  template.vl <- seq(100, 200, length.out = 10)
  template.vl <- sort(unique(round(template.vl)))

  res <- scaleTemplate(template, template.vl)

  dres.actual <- unlist(lapply(res, function(res.i) unique(lapply(res.i, length))))
  dres.exp    <- template.vl

  expect_equal(dres.actual, dres.exp)
})


test_that("Testing correctness of ordering of elements in output list", {

  set.seed(20180927)

  template    <- list(rnorm(100), rnorm(110), rnorm(120), rnorm(130), rnorm(140))
  template.vl <- seq(100, 200, length.out = 10)
  template.vl <- sort(unique(round(template.vl)))

  res <- scaleTemplate(template, template.vl)

  dres.actual <- unlist(lapply(res, function(res.i) res.i[[1]][1]))
  dres.exp <- c(-0.71691806384448, -0.894959652480998, -0.893147165197797,
                -0.899933244566838, -0.89138096343759, -0.900057124878235, -0.902450851159602,
                -0.898388774448527, -0.895921300251811, -0.892782809450523)
  expect_equal(dres.actual, dres.exp)

})

