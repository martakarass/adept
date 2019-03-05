
# rm(list = ls())

context("Testing segmentPattern function.")

######################
######################
######################
###   Examples 1   ###
######################
######################
######################

test_that("Example 1(a): no noise in signal, all pattern occurences of the same
          length (vector length: 101)", {

  ## Generate signal and template
  x0 <- cos(seq(0, 2 * pi * 10, length.out = 1001))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 110)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- sort(out$tau_i)
  res.exp <- c(1, 101, 201, 302, 401, 501, 601, 701, 801, 901)
  expect_equal(res, res.exp)

  res <- sort(out$T_i)
  res.exp <- c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100)
  expect_equal(res, res.exp)

  res <- sort(out$sim_i)
  res.exp <- c(0.99947934173182, 0.999479341731821, 0.999479341731821, 0.999479341731821,
               0.999479341731821, 0.999479341731821, 0.999479341731821, 0.999479341731821,
               0.999479341731821, 0.999479341731821)
  expect_equal(res, res.exp)

})



test_that("Example 1(b): no noise in signal, all pattern occurences of the same length
          (vector length: 101), include in grid of considered pattern durations the
          true value 101", {

  ## Generate signal and template
  x0 <- cos(seq(0, 2 * pi * 10, length.out = 1001))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 101, 102, 110)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- sort(out$tau_i)
  res.exp <- c(1, 101, 201, 301, 401, 501, 601, 701, 801, 901)
  expect_equal(res, res.exp)

  res <- sort(out$T_i)
  res.exp <- c(101, 101, 101, 101, 101, 101, 101, 101, 101, 101)
  expect_equal(res, res.exp)

  res <- sort(out$sim_i)
  res.exp <- c(0.999999999999999, 0.999999999999999, 1, 1, 1, 1, 1, 1, 1,
               1)
  expect_equal(res, res.exp)
})



test_that("Example 1(c): no noise in signal, all pattern occurences of the same length
          (vector length: 101), include in grid of considered pattern durations the
          true value 101,
          modify assumed signal frequency and scale of grid of pattern durations", {

  ## Generate signal and template
  x0 <- cos(seq(0, 2 * pi * 10, length.out = 1001))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 101, 102, 110) * 0.01
  out <- segmentPattern(x = x,
                        x.fs = 100,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- sort(out$tau_i)
  res.exp <- c(1, 101, 201, 301, 401, 501, 601, 701, 801, 901)
  expect_equal(res, res.exp)

  res <- sort(out$T_i)
  res.exp <- c(101, 101, 101, 101, 101, 101, 101, 101, 101, 101)
  expect_equal(res, res.exp)

  res <- sort(out$sim_i)
  res.exp <- c(0.999999999999999, 0.999999999999999, 1, 1, 1, 1, 1, 1, 1,
               1)
  expect_equal(res, res.exp)
})



######################
######################
######################
###   Examples 2   ###
######################
######################
######################

# set.seed(1)
# s.grid <- sample(60:120, size = 10)
## changed to correct for sample function bias in new release of R
s.grid <- c(76L, 82L, 93L, 112L, 71L, 110L, 111L, 95L, 118L, 63L)

if (exists("template")) rm(template)

test_that("Example 2(a): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length", {

  ## Grid of different true pattern occurence durations
  template <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }

  ## Assume grid of assumed duration of pattern occurrences
  pattern.dur.seq <- 60:120
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- out$tau_i
  res.exp <- c(1, 76, 157, 249, 360, 430, 539, 649, 743, 860)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(76, 82, 93, 112, 71, 110, 111, 95, 118, 63)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.999999999999989, 0.999999999999991, 0.99999999999999, 0.99999999999999,
               0.999999999999991, 0.999999999999988, 0.999999999999991, 0.999999999999991,
               0.999999999999989, 0.999999999999992)
  expect_equal(res, res.exp)
})




test_that("Example 2(b): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; What happens if we define very
          'poor' grid of assumed duration of pattern occurrences?", {

  ## Grid of different true pattern occurence durations
  template <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }

  ## Assume grid of assumed duration of pattern occurrences
  pattern.dur.seq <- c(60, 80, 100, 120)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- out$tau_i
  res.exp <- c(9, 77, 163, 245, 366, 435, 534, 656, 742, 861)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(60, 80, 80, 120, 60, 100, 120, 80, 120, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.991353694993282, 0.999831828958601, 0.994900605403087, 0.998088931040552,
               0.993303282206735, 0.998018011681711, 0.997327871634614, 0.993913153322669,
               0.999911747300496, 0.997997165500962)
  expect_equal(res, res.exp)
})




# test_that("Example 2(c): no noise in signal (no peak fine-tuning employed),
#           pattern occurences of different length; Use covariance, and even more
#           'poor' grid of assumed duration of pattern occurrences", {
#
#   ## Grid of different true pattern occurence durations
#   template <- cos(seq(0, 2 * pi, length.out = 200))
#
#   ## Generate signal x that consists of "glued" pattern occurrences of different length
#   x <- numeric()
#   for (ss in s.grid){
#     templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
#     if (length(x)>0){
#       x <- x[-length(x)]
#     }
#     x <- c(x, templ0)
#   }
#
#   ## Assume grid of assumed duration of pattern occurrences
#   pattern.dur.seq <- c(60, 120)
#   out <- segmentPattern(x = x,
#                         x.fs = 1,
#                         template = template,
#                         pattern.dur.seq = pattern.dur.seq,
#                         similarity.measure = "cov")
#
#   res <- out$tau_i
#   res.exp <- c(9, 87, 174, 245, 366, 455, 534, 666, 742, 861)
#   expect_equal(res, res.exp)
#
#   res <- out$T_i
#   res.exp <- c(60, 60, 60, 120, 60, 60, 120, 60, 120, 60)
#   expect_equal(res, res.exp)
#
#   res <- out$sim_i
#   res.exp <- c(0.591863061841216, 0.540762733635577, 0.455917527932434, 0.728925513944514,
#                0.634801827775348, 0.35289880518215, 0.730662908614629, 0.442073525854687,
#                0.718248707449562, 0.69894752810577)
#   expect_equal(res, res.exp)
# })



test_that("Example 2(c): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; Use covariance", {

  ## Grid of different true pattern occurence durations
  template <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }

  ## Assume grid of assumed duration of pattern occurrences
  pattern.dur.seq <- c(60, 80, 100, 120)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cov")

  res <- out$tau_i
  res.exp <- c(1, 87, 154, 255, 355, 435, 534, 656, 742, 861)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(80, 60, 100, 100, 80, 100, 120, 80, 120, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.719501786387912, 0.540762733635577, 0.731177684709071, 0.662370934346171,
               0.737357140142741, 0.672095166679452, 0.730662908614629, 0.629670599324224,
               0.718248707449562, 0.69894752810577)
  expect_equal(res, res.exp)
  })




######################
######################
######################
###   Examples 3   ###
######################
######################
######################



test_that("Example 3(a): add noise in signal generation,
          pattern occurences of different length in singal generation", {

  template <- cos(seq(0, 2 * pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }

  set.seed(1)
  x <- x + rnorm(length(x), sd = 0.3)

  pattern.dur.seq <- seq(60, 100, by = 5)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- out$tau_i
  res.exp <- c(4, 73, 162, 255, 362, 435, 545, 645, 752, 862)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(70, 90, 85, 100, 65, 100, 100, 100, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.932841298514519, 0.937485762099541, 0.903795488242938, 0.902713349912544,
               0.915232790752724, 0.87794212230542, 0.894118326687643, 0.907734430293421,
               0.893770550717209, 0.93223222750356)
  expect_equal(res, res.exp)
})




test_that("Example 3(b): add noise in signal generation,
          pattern occurences of different length in singal generation", {

  ## Generate signal and template
  template <- cos(seq(0, 2 * pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }

  set.seed(1)
  x <- x + rnorm(length(x), sd = 0.3)

  ## Without maxima detection
  pattern.dur.seq <- seq(60, 100, by = 5)


  # similarity.measure.thresh = 0.0
  # x.adept.ma.W = NULL
  # finetune = NULL
  # finetune.maxima.ma.W = NULL
  # finetune.maxima.nbh.W = NULL
  # run.parallel = FALSE
  # run.parallel.cores = NULL
  # x.cut = TRUE
  # x.cut.vl = 6000
  # compute.template.idx = FALSE

  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- out$tau_i
  res.exp <- c(4, 73, 162, 255, 362, 435, 545, 645, 752, 862)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(70, 90, 85, 100, 65, 100, 100, 100, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.932841298514519, 0.937485762099541, 0.903795488242938, 0.902713349912544,
               0.915232790752724, 0.87794212230542, 0.894118326687643, 0.907734430293421,
               0.893770550717209, 0.93223222750356)
  expect_equal(res, res.exp)


  # similarity.measure = "cor"
  # finetune = "maxima"
  # finetune.maxima.ma.W = 30
  # finetune.maxima.nbh.W = 120

  ## Use maxima detection (for comparison)
  pattern.dur.seq <- seq(60, 100, by = 5)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        finetune = "maxima",
                        finetune.maxima.ma.W = 30,
                        finetune.maxima.nbh.W = 120)

  res <- out$tau_i
  res.exp <- c(15, 75, 160, 253, 355, 435, 542, 648, 750, 849)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(61, 86, 92, 100, 77, 100, 100, 98, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.934034172881348, 0.937485762099541, 0.914074642973907, 0.902713349912544,
               0.915232790752724, 0.87794212230542, 0.894118326687643, 0.907734430293421,
               0.888048574952028, 0.93223222750356)
  expect_equal(res, res.exp)



  ## Use maxima detection, use thresholding on covariance
  pattern.dur.seq <- seq(60, 100, by = 5)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure.thresh = 0.7,
                        similarity.measure = "cov",
                        finetune = "maxima",
                        finetune.maxima.ma.W = 30,
                        finetune.maxima.nbh.W = 120)

  res <- out$tau_i
  res.exp <- c(15, 160, 648, 849)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(61, 92, 98, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.724850738169134, 0.736528364289848, 0.717833003791632, 0.743334116677552)
  expect_equal(res, res.exp)



  ## Use maxima detection, use thresholding on covariance,
  ## and add smoothing
  pattern.dur.seq <- seq(60, 100, by = 5)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure.thresh = -1,
                        x.adept.ma.W = 30,
                        similarity.measure = "cov",
                        finetune = "maxima",
                        finetune.maxima.ma.W = 30,
                        finetune.maxima.nbh.W = 120)

  res <- out$tau_i
  res.exp <- c(15, 75, 160, 253, 355, 435, 542, 648, 751, 857)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(61, 86, 92, 100, 77, 100, 100, 98, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.398855621994664, 0.574913156763604, 0.613893094434027, 0.562563743167577,
               0.544951870561808, 0.562680022332807, 0.563855908578859, 0.614970735584006,
               0.593536450844769, 0.320552985778579)
  expect_equal(res, res.exp)


  ## increase range of allowed pattern duration (pattern.dur.seq)
  pattern.dur.seq <- seq(50, 150, by = 5)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        finetune = "maxima",
                        finetune.maxima.ma.W = 30,
                        finetune.maxima.nbh.W = 120)

  res <- out$tau_i
  res.exp <- c(15, 75, 160, 251, 355, 431, 541, 648, 745, 857)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(61, 86, 92, 105, 77, 111, 108, 98, 113, 52)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.934034172881348, 0.937485762099541, 0.914074642973907, 0.920834847349191,
               0.915232790752724, 0.896116452119818, 0.91682323488697, 0.9072682178912,
               0.898872120860981, 0.93223222750356)
  expect_equal(res, res.exp)

})






######################
######################
######################
###   Examples 4   ###
######################
######################
######################

## CHECKING FOR x.cut = TRUE correctness


test_that("Example 4(a): no noise in signal, all pattern occurences of the same
          length, long signal. We check for same behaviour with and without x.cut = TRUE", {

  ## Generate signal and template
  ## approx 30 min for data collected at freq
  x0 <- cos(seq(0, 2 * pi * 10 * 20, length.out = 1000 * 20 + 1))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 101)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        x.cut = TRUE,
                        x.cut.vl = 6000)

  res <- as.vector(summary(out$tau_i))
  res.exp <- c(1, 4976, 9951, 9951, 14926, 19901)
  expect_equal(res, res.exp)

  res <- as.vector(summary(out$T_i))
  res.exp <- c(101, 101, 101, 101, 101, 101)
  expect_equal(res, res.exp)

  res <- as.vector(summary(out$sim_i))
  res.exp <- c(0.999999999999998, 1, 1, 1, 1, 1)
  expect_equal(res, res.exp)

  ## compare with no cut
  out2 <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        x.cut = FALSE)

  expect_equal(out$tau_i, out2$tau_i)
  expect_equal(out$T_i,   out2$T_i)
  expect_equal(out$sim_i, out2$sim_i)
})



## TODO

test_that("Example 4(b): no noise in signal, all pattern occurences of the same
          length, long signal. We check for same behaviour with and without x.cut = TRUE", {

  ## Grid of different true pattern occurence durations
  template1 <- cos(seq(0, 2 * pi, length.out = 200))
  template2 <- c(rev(seq(-1, 1, length.out = 100)), seq(-1, 1, length.out = 100))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in as.vector(replicate(1, seq(80, 120)))){
    ## Add piece from template1
    templ0 <- approx(seq(0, 1, length.out = 200), template1, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0) x <- x[-length(x)]
    x <- c(x, templ0)
    ## Add piece from template2
    templ0 <- approx(seq(0, 1, length.out = 200), template2, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0) x <- x[-length(x)]
    x <- c(x, templ0)
  }

  pattern.dur.seq <- seq(80, 120, by = 4)

  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = list(template1, template2),
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cov",
                        x.cut = TRUE,
                        x.cut.vl = 6000,
                        compute.template.idx = TRUE)

  out2 <- segmentPattern(x = x,
                        x.fs = 1,
                        template = list(template1, template2),
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cov",
                        x.cut = FALSE,
                        compute.template.idx = TRUE)

  expect_equal(out$tau_i, out2$tau_i)
  expect_equal(out$T_i,   out2$T_i)
  expect_equal(out$sim_i, out2$sim_i)
  expect_equal(out$template_i, out2$template_i)
})


# data.frame(out$tau_i, out2$tau_i, out$T_i, out2$T_i, out$sim_i, out2$sim_i)



######################
######################
######################
###   Examples 5   ###
######################
######################
######################

## PARALLEL COMPUTATION CHECKING


test_that("Example 5(a). Setting: no noise in signal, all pattern occurences of the same
          length, long signal. We check for the same behaviour with and without parallel
          computation with fixed number of 2 workers.", {

  ## Generate signal and template
  ## approx 30 min for data collected at freq
  x0 <- cos(seq(0, 2 * pi * 10 * 20, length.out = 1000 * 20 + 1))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 101)

  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        x.cut = TRUE,
                        x.cut.vl = 6000,
                        compute.template.idx = TRUE)

  out2 <- segmentPattern(x = x,
                         x.fs = 1,
                         template = template,
                         pattern.dur.seq = pattern.dur.seq,
                         similarity.measure = "cor",
                         x.cut = TRUE,
                         x.cut.vl = 6000,
                         run.parallel = TRUE,
                         run.parallel.cores = 2,
                         compute.template.idx = TRUE)

  expect_equal(out$tau_i, out2$tau_i)
  expect_equal(out$T_i,   out2$T_i)
  expect_equal(out$sim_i, out2$sim_i)
  expect_equal(out$template_i, out2$template_i)
})



test_that("Example 5(b). Setting: no noise in signal, all pattern occurences of the same
          length, long signal. We check for the same behaviour with and without parallel
          computation with maximum number of workers available.", {

  ## Generate signal and template
  ## approx 30 min for data collected at freq
  x0 <- cos(seq(0, 2 * pi * 10 * 20, length.out = 1000 * 20 + 1))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 101)

  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        x.cut = TRUE,
                        x.cut.vl = 6000,
                        compute.template.idx = TRUE)

  out2 <- segmentPattern(x = x,
                         x.fs = 1,
                         template = template,
                         pattern.dur.seq = pattern.dur.seq,
                         similarity.measure = "cor",
                         x.cut = TRUE,
                         x.cut.vl = 6000,
                         run.parallel = TRUE,
                         run.parallel.cores = NULL,
                         compute.template.idx = TRUE)

  expect_equal(out$tau_i, out2$tau_i)
  expect_equal(out$T_i,   out2$T_i)
  expect_equal(out$sim_i, out2$sim_i)
  expect_equal(out$template_i, out2$template_i)
})



test_that("Example 5(c). Setting: no noise in signal, all pattern occurences of the same
          length, long signal. We check for the same behaviour with and without parallel
          computation with (maximum-1) number of workers available.", {

  ## Generate signal and template
  ## approx 30 min for data collected at freq
  x0 <- cos(seq(0, 2 * pi * 10 * 20, length.out = 1000 * 20 + 1))
  x  <- x0
  template <- x0[1:101]

  pattern.dur.seq <- c(90, 100, 101)

  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor",
                        x.cut = TRUE,
                        x.cut.vl = 6000,
                        compute.template.idx = TRUE)

  out2 <- segmentPattern(x = x,
                         x.fs = 1,
                         template = template,
                         pattern.dur.seq = pattern.dur.seq,
                         similarity.measure = "cor",
                         x.cut = TRUE,
                         x.cut.vl = 6000,
                         run.parallel = TRUE,
                         run.parallel.cores = availableCores(),
                         compute.template.idx = TRUE)

  expect_equal(out$tau_i, out2$tau_i)
  expect_equal(out$T_i,   out2$T_i)
  expect_equal(out$sim_i, out2$sim_i)
  expect_equal(out$template_i, out2$template_i)
})




######################
######################
######################
###   Examples 6   ###
######################
######################
######################

## CHECKING FOR TEMPLATE INDEX

# set.seed(1)
# s.grid2 <- sample(60:120, size = 50)
s.grid2 <- c(76L, 82L, 93L, 112L, 71L, 110L, 111L, 95L, 118L, 63L, 70L,
            68L, 117L, 78L, 96L, 119L, 92L, 103L, 120L, 104L, 98L, 109L,
            85L, 64L, 69L, 73L, 60L, 113L, 88L, 115L, 74L, 77L, 90L, 65L,
            105L, 89L, 79L, 62L, 102L, 106L, 100L, 72L, 101L, 81L, 80L, 84L,
            94L, 66L, 86L, 99L)

test_that("Example 6(a): test returning template index matrix", {


  ## Grid of different true pattern occurence durations
  template1 <- cos(seq(0, 2 * pi, length.out = 200))
  template2 <- c(rev(seq(-1, 1, length.out = 100)), seq(-1, 1, length.out = 100))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid2){
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
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = list(template1, template2),
                        pattern.dur.seq = sort(s.grid2),
                        similarity.measure = "cor",
                        compute.template.idx = TRUE)

  expect_true(all(out$template_i == rep(c(1,2), 50)))
})






