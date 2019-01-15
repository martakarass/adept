

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

if (exists("template")) rm(template)

test_that("Example 2(a): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
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
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
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




test_that("Example 2(c): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; Use covariance, and even more
          'poor' grid of assumed duration of pattern occurrences", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
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
  pattern.dur.seq <- c(60, 120)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cov")

  res <- out$tau_i
  res.exp <- c(9, 87, 174, 245, 366, 455, 534, 666, 742, 861)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(60, 60, 60, 120, 60, 60, 120, 60, 120, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.591863061841216, 0.540762733635577, 0.455917527932434, 0.728925513944514,
               0.634801827775348, 0.35289880518215, 0.730662908614629, 0.442073525854687,
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

  ## Generate signal and template
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  template <- cos(seq(0, 2 * pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }
  x <- x + rnorm(length(x), sd = 0.3)

  pattern.dur.seq <- seq(60, 100, by = 5)
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- out$tau_i
  res.exp <- c(4, 73, 162, 254, 358, 435, 545, 646, 752, 862)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(70, 90, 85, 100, 75, 100, 100, 100, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.926292653865246, 0.943177988579057, 0.909149333113241, 0.905514398680946,
               0.922521494082868, 0.878835859702908, 0.896399820840297, 0.908359071401425,
               0.892506179992972, 0.932539440191918)
  expect_equal(res, res.exp)
})




test_that("Example 3(b): add noise in signal generation,
          pattern occurences of different length in singal generation", {

  ## Generate signal and template
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  template <- cos(seq(0, 2 * pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), template, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }
  x <- x + rnorm(length(x), sd = 0.3)



  ## Without maxima detection
  pattern.dur.seq <- seq(60, 100, by = 5)


  similarity.measure.thresh = 0.0
  x.adept.ma.W = NULL
  finetune = NULL
  finetune.maxima.ma.W = NULL
  finetune.maxima.nbh.W = NULL
  run.parallel = FALSE
  run.parallel.cores = NULL
  x.cut = TRUE
  x.cut.vl = 6000
  compute.template.idx = FALSE

  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = template,
                        pattern.dur.seq = pattern.dur.seq,
                        similarity.measure = "cor")

  res <- out$tau_i
  res.exp <- c(4, 73, 162, 254, 358, 435, 545, 646, 752, 862)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(70, 90, 85, 100, 75, 100, 100, 100, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.926292653865246, 0.943177988579057, 0.909149333113241, 0.905514398680946,
               0.922521494082868, 0.878835859702908, 0.896399820840297, 0.908359071401425,
               0.892506179992972, 0.932539440191918)
  expect_equal(res, res.exp)


  similarity.measure = "cor"
  finetune = "maxima"
  finetune.maxima.ma.W = 30
  finetune.maxima.nbh.W = 120

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
  res.exp <- c(15, 77, 160, 257, 356, 434, 544, 651, 750, 849)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(63, 84, 89, 100, 79, 100, 100, 97, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.934011470294625, 0.943177988579057, 0.917625765654425, 0.905514398680946,
               0.922521494082868, 0.878835859702908, 0.896399820840297, 0.908359071401425,
               0.89029293486179, 0.932539440191918)
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
  res.exp <- c(15, 77, 160, 651, 849)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(63, 84, 89, 97, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.705380000974917, 0.735921165890131, 0.759885566360049, 0.727241676572737,
               0.74274381619032)
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
  res.exp <- c(15, 77, 160, 257, 356, 434, 544, 651, 751, 852)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(63, 84, 89, 100, 79, 100, 100, 97, 100, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.387428828230924, 0.593344449562667, 0.628648114234425, 0.564177922975735,
               0.524566992985259, 0.57315100219008, 0.561412683480699, 0.61036916804888,
               0.589173737378179, 0.334041868743848)
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
  res.exp <- c(15, 77, 160, 248, 356, 434, 537, 651, 747, 852)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(63, 84, 89, 109, 79, 104, 115, 97, 106, 57)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.934011470294625, 0.943177988579057, 0.917625765654425, 0.905514398680946,
               0.922521494082868, 0.878835859702908, 0.913581497874648, 0.902543119926644,
               0.892506179992972, 0.932539440191918)
  expect_equal(res, res.exp)

})






######################
######################
######################
###   Examples 4   ###
######################
######################
######################



test_that("Example 4(a): no noise in signal, all pattern occurences of the same
          length, LONG SIGNAL", {

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

  ## compare with cut, parallel
  out3 <- segmentPattern(x = x,
                         x.fs = 1,
                         template = template,
                         pattern.dur.seq = pattern.dur.seq,
                         similarity.measure = "cor",
                         x.cut = TRUE,
                         x.cut.vl = 6000,
                         run.parallel = TRUE,
                         run.parallel.cores = 2)


  expect_equal(out$tau_i, out2$tau_i)
  expect_equal(out$tau_i, out3$tau_i)

  expect_equal(out$T_i, out2$T_i)
  expect_equal(out$T_i, out3$T_i)

  expect_equal(out$sim_i, out2$sim_i)
  expect_equal(out$sim_i, out3$sim_i)
})




######################
######################
######################
###   Examples 4   ###
######################
######################
######################


test_that("Example 5(a): test returning template index matrix", {


  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 50)
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
  out <- segmentPattern(x = x,
                        x.fs = 1,
                        template = list(template1, template2),
                        pattern.dur.seq = sort(s.grid),
                        similarity.measure = "cor",
                        compute.template.idx = TRUE)

  expect_true(all(out$template_i == rep(c(1,2), 50)))
})






