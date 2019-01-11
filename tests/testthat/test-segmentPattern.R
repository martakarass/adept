

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


test_that("Example 2(a): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  templ <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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

  res <- sort(out$tau_i)
  res.exp <- c(1, 76, 157, 249, 360, 430, 539, 649, 743, 860)
  expect_equal(res, res.exp)

  res <- sort(out$T_i)
  res.exp <- c(63, 71, 76, 82, 93, 95, 110, 111, 112, 118)
  expect_equal(res, res.exp)

  res <- sort(out$sim_i)
  res.exp <- c(0.999999978863449, 0.999999986147232, 0.999999987426459, 0.999999987616627,
               0.999999987957358, 0.999999988038671, 0.999999988177003, 0.999999988262704,
               0.999999988270857, 0.999999988274299)
  expect_equal(res, res.exp)
})




test_that("Example 2(b): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; What happens if we define very
          'poor' grid of assumed duration of pattern occurrences?", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  templ <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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
  res.exp <- c(9, 77, 163, 245, 366, 435, 534, 657, 742, 861)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(60, 80, 80, 120, 60, 100, 120, 80, 120, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.991356582934436, 0.999831987402593, 0.994901487934483, 0.998088156917189,
              0.993305509171119, 0.998023479621162, 0.997327183352651, 0.993914137573433,
              0.999911593085322, 0.99799791084828)
  expect_equal(res, res.exp)
})




test_that("Example 2(c): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; Use covariance, and even more
          'poor' grid of assumed duration of pattern occurrences", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  templ <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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
  res.exp <- c(9, 87, 174, 245, 365, 455, 534, 666, 742, 862)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(60, 60, 60, 120, 60, 60, 120, 60, 120, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.59186478601464, 0.540764650172252, 0.455919521310586, 0.728924948585804,
               0.634803250987258, 0.352900628868909, 0.730662404365108, 0.442075511219409,
               0.718248596674937, 0.698948050109923)
  expect_equal(res, res.exp)
})





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


test_that("Example 2(a): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length", {

            ## Grid of different true pattern occurence durations
            set.seed(1)
            s.grid <- sample(60:120, size = 10)
            templ <- cos(seq(0, 2 * pi, length.out = 200))

            ## Generate signal x that consists of "glued" pattern occurrences of different length
            x <- numeric()
            for (ss in s.grid){
              templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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

            res <- sort(out$tau_i)
            res.exp <- c(1, 76, 157, 249, 360, 430, 539, 649, 743, 860)
            expect_equal(res, res.exp)

            res <- sort(out$T_i)
            res.exp <- c(63, 71, 76, 82, 93, 95, 110, 111, 112, 118)
            expect_equal(res, res.exp)

            res <- sort(out$sim_i)
            res.exp <- c(0.999999978863449, 0.999999986147232, 0.999999987426459, 0.999999987616627,
                         0.999999987957358, 0.999999988038671, 0.999999988177003, 0.999999988262704,
                         0.999999988270857, 0.999999988274299)
            expect_equal(res, res.exp)
            })




test_that("Example 2(b): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; What happens if we define very
          'poor' grid of assumed duration of pattern occurrences?", {

            ## Grid of different true pattern occurence durations
            set.seed(1)
            s.grid <- sample(60:120, size = 10)
            templ <- cos(seq(0, 2 * pi, length.out = 200))

            ## Generate signal x that consists of "glued" pattern occurrences of different length
            x <- numeric()
            for (ss in s.grid){
              templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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
            res.exp <- c(9, 77, 163, 245, 366, 435, 534, 657, 742, 861)
            expect_equal(res, res.exp)

            res <- out$T_i
            res.exp <- c(60, 80, 80, 120, 60, 100, 120, 80, 120, 60)
            expect_equal(res, res.exp)

            res <- out$sim_i
            res.exp <- c(0.991356582934436, 0.999831987402593, 0.994901487934483, 0.998088156917189,
                         0.993305509171119, 0.998023479621162, 0.997327183352651, 0.993914137573433,
                         0.999911593085322, 0.99799791084828)
            expect_equal(res, res.exp)
            })




test_that("Example 2(d): no noise in signal (no peak fine-tuning employed),
          pattern occurences of different length; Use covariance, and even more
          'poor' grid of assumed duration of pattern occurrences;
          Employ threshold for covariance", {

  ## Grid of different true pattern occurence durations
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  templ <- cos(seq(0, 2 * pi, length.out = 200))

  ## Generate signal x that consists of "glued" pattern occurrences of different length
  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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
                        similarity.measure = "cov",
                        similarity.measure.thresh = 0.6)

  res <- out$tau_i
  res.exp <- c(245, 365, 534, 742, 862)
  expect_equal(res, res.exp)

  res <- out$T_i
  res.exp <- c(120, 60, 120, 120, 60)
  expect_equal(res, res.exp)

  res <- out$sim_i
  res.exp <- c(0.728924948585804, 0.634803250987258, 0.730662404365108, 0.718248596674937,
               0.698948050109923)
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
  templ <- cos(seq(0, 2 * pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
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
  res.exp <- c(0.926291306033366, 0.943178533162363, 0.909156293953748, 0.905525014013764,
               0.922523794573569, 0.878846405199884, 0.896401753440613, 0.908357873425458,
               0.892508399360682, 0.932536508269999)
  expect_equal(res, res.exp)
})




test_that("Example 3(b): add noise in signal generation,
          pattern occurences of different length in singal generation", {

  ## Generate signal and template
  set.seed(1)
  s.grid <- sample(60:120, size = 10)
  templ <- cos(seq(0, 2 * pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }
  x <- x + rnorm(length(x), sd = 0.3)



  ## Without maxima detection
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
  res.exp <- c(0.926291306033366, 0.943178533162363, 0.909156293953748, 0.905525014013764,
               0.922523794573569, 0.878846405199884, 0.896401753440613, 0.908357873425458,
               0.892508399360682, 0.932536508269999)
  expect_equal(res, res.exp)



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
  res.exp <- c(0.934008825205221, 0.943178533162363, 0.917627061535885, 0.905525014013764,
               0.922523794573569, 0.878846405199884, 0.896401753440613, 0.908357873425458,
               0.890287715999962, 0.932536508269999)
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
  res.exp <- c(0.705378003362227, 0.735918052569151, 0.759884581289357, 0.727240717460776,
               0.742750320506588)
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
  res.exp <- c(0.387425918617852, 0.59334391287926, 0.628645307416092, 0.564182215998019,
               0.524567526202169, 0.573154721743202, 0.561415062231443, 0.610366994688926,
               0.589176030015489, 0.33404216471104)
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
  res.exp <- c(0.934008825205221, 0.943178533162363, 0.917627061535885, 0.905525014013764,
               0.922523794573569, 0.878846405199884, 0.913586332840362, 0.902548497003548,
               0.892508399360682, 0.932536508269999)
  expect_equal(res, res.exp)

})



