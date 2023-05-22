

context("Testing maxAndTune.")

test_that("Sute 1. Test if maxAndTune works correctly. Case: no noise in signal, no tunning, one of the pattern durations is perfectly.", {

  ## Params
  x0 <- cos(seq(0, 2 * pi * 10, length.out = 1001))
  x  <- x0
  template <- list(x0[1:101])
  pattern.dur.seq <- c(80, 90, 100)

  finetune <- NULL
  finetune.maxima.x <- NULL
  finetune.maxima.nbh.vl <- NULL
  similarity.measure.thresh <- 0.0

  ## Partial objects
  x.fs <- 1
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)
  x.smoothed <- x
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cov")$similarity

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  # expect_equal(sort(out$tau_i),
  #              c(1, 101, 201, 302, 401, 501, 601, 701, 801, 901))
  # expect_equal(sort(out$T_i),
  #              c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100))
  ## @MK May 11, 2019: account for micro differences in similarity matrix
  ## computation osx/win which may yield up to 1 index difference in stride location
  ## estimation
  expect_true(max(abs(sort(out$tau_i) - c(1, 101, 201, 302, 401, 501, 601, 701, 801, 901))) <= 1)
  expect_true(max(abs(sort(out$T_i) - c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100))) <= 1)
})


test_that("Sute 2. Test if maxAndTune works correctly. Case: no noise in signal, no tunning, one of the pattern durations is perfectly.", {

  ## Params
  x0 <- cos(seq(0, 2 * pi * 10, length.out = 1001))
  x  <- x0
  template <- list(x0[1:90])
  pattern.dur.seq <- c(80, 90, 100)

  finetune <- NULL
  finetune.maxima.x <- NULL
  finetune.maxima.nbh.vl <- NULL
  similarity.measure.thresh <- 0.0

  ## Partial objects
  x.fs <- 1
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)
  x.smoothed <- x
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")$similarity

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_equal(sort(out$tau_i),
               c(1, 101, 201, 301, 401, 501, 601, 701, 801, 901))
  expect_equal(sort(out$T_i),
               c(90, 90, 90, 90, 90, 90, 90, 90, 90, 90))

})



## Replace to avoid sample() function bias
# s.grid.sample <- sample(c(80, 100, 120), size = 50, replace = TRUE)
s.grid.sample <- c(120, 100, 80, 120, 120, 80, 120, 80, 100, 80, 100, 120, 120,
                   100, 120, 100, 80, 120, 80, 100, 80, 100, 80, 120, 120, 120,
                   120, 100, 100, 80, 120, 100, 120, 80, 80, 80, 120, 100, 100,
                   100, 100, 80, 80, 120, 120, 80, 100, 120, 80, 100)

test_that("Sute 3. Test if maxAndTune works correctly. Case: no noise in signal, no tunning, one of the pattern durations is perfectly.", {

  set.seed(20181019)
  templ <- sin(seq(0, pi, length.out = 200))

  x <- numeric()
  # for (ss in sample(s.grid, size = 50, replace = TRUE)){
  for (ss in s.grid.sample){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }
  # x <- x + rnorm(length(x), sd = 0.4)

  template <- list(templ)
  pattern.dur.seq <- c(80, 100, 120)

  finetune <- NULL
  finetune.maxima.x <- NULL
  finetune.maxima.nbh.vl <- NULL
  similarity.measure.thresh <- 0.0

  ## Partial objects
  x.fs <- 1
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)
  x.smoothed <- x
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")$similarity

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_true(all(sort(out$tau_i) == sort(which(x == 0))))

  expect_equal(sort(out$tau_i),
               c(1, 120, 219, 298, 417, 536, 615, 734, 813, 912, 991, 1090,
                 1209, 1328, 1427, 1546, 1645, 1724, 1843, 1922, 2021, 2100, 2199,
                 2278, 2397, 2516, 2635, 2754, 2853, 2952, 3031, 3150, 3249, 3368,
                 3447, 3526, 3605, 3724, 3823, 3922, 4021, 4120, 4199, 4278, 4397,
                 4516, 4595, 4694, 4813, 4892))
  expect_equal(sort(out$T_i),
               c(80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80,
                 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                 100, 100, 100, 100, 120, 120, 120, 120, 120, 120, 120, 120, 120,
                 120, 120, 120, 120, 120, 120, 120, 120, 120))
  expect_equal(sort(out$sim_i),
               c(0.999999999999999, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))

})





test_that("Sute 4. Test if maxAndTune works correctly. Case: add noise to signal, no tunning", {

  set.seed(20181019)
  templ <- sin(seq(0, pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid.sample){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }
  x <- x + rnorm(length(x), sd = 0.4)

  template <- list(templ)
  pattern.dur.seq <- c(70, 80, 100, 120)

  finetune <- NULL
  finetune.maxima.x <- NULL
  finetune.maxima.nbh.vl <- NULL
  similarity.measure.thresh <- 0.0

  ## Partial objects
  x.fs <- 1
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)
  x.smoothed <- x
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")$similarity

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_equal(sort(out$tau_i),
               c(2, 121, 226, 298, 423, 522, 635, 736, 815, 916, 988, 1087,
                 1208, 1327, 1429, 1548, 1647, 1735, 1840, 1919, 2025, 2110, 2202,
                 2281, 2392, 2519, 2638, 2771, 2852, 2955, 3038, 3137, 3256, 3443,
                 3530, 3601, 3724, 3826, 3945, 4028, 4120, 4211, 4280, 4401, 4584,
                 4703, 4814, 4912))
  expect_equal(sort(out$T_i),
               c(70, 70, 70, 70, 70, 70, 80, 80, 80, 80, 80, 80, 80, 80, 80,
                 80, 80, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                 100, 100, 100, 100, 100, 100, 120, 120, 120, 120, 120, 120, 120,
                 120, 120, 120, 120, 120, 120, 120))
  expect_equal(sort(out$sim_i),
               c(0.408537323284492, 0.435800410963697, 0.452862668704024, 0.456536259437437,
                 0.470256190840373, 0.473409842060613, 0.47914175824858, 0.48000754936811,
                 0.528332523105161, 0.541308746577531, 0.545290195435071, 0.550199791127224,
                 0.562384888929942, 0.569030435117554, 0.57375553093907, 0.577748913120457,
                 0.578898693945664, 0.598093686855695, 0.604150050986692, 0.6052739242771,
                 0.610176794203946, 0.612023739225788, 0.614492285215261, 0.620227553733224,
                 0.621164966971297, 0.622783036449866, 0.627319873312495, 0.633248513173231,
                 0.635107424444073, 0.635521689774035, 0.636802853348891, 0.647105239817411,
                 0.653594865973984, 0.65658700595948, 0.661407919988764, 0.66171351823905,
                 0.669422058643667, 0.677205987228789, 0.677569883368504, 0.681830754313774,
                 0.682351757248779, 0.682661399905778, 0.685176833969805, 0.698700404217461,
                 0.70183423417834, 0.743425396253583, 0.755487202273469, 0.764657255410779
               ))

})



test_that("Sute 5. Test if maxAndTune works correctly. Case: add noise to signal, no tunning", {

  set.seed(20181019)
  templ <- sin(seq(0, pi, length.out = 200))

  x <- numeric()
  for (ss in s.grid.sample){
    templ0 <- approx(seq(0, 1, length.out = 200), templ, xout = seq(0, 1, length.out = ss))$y
    if (length(x)>0){
      x <- x[-length(x)]
    }
    x <- c(x, templ0)
  }
  x <- x + rnorm(length(x), sd = 0.4)

  template <- list(templ)
  pattern.dur.seq <- c(70, 80, 100, 120)

  finetune <- NULL
  finetune.maxima.x <- NULL
  finetune.maxima.nbh.vl <- NULL
  similarity.measure.thresh <- 0.6

  ## Partial objects
  x.fs <- 1
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)
  x.smoothed <- x
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")$similarity

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_equal(sort(out$tau_i),
               c(2, 121, 298, 522, 635, 736, 815, 1087, 1327, 1429, 1548, 1647,
                 1735, 1840, 1919, 2025, 2202, 2392, 2519, 2638, 2852, 2955, 3137,
                 3256, 3443, 3826, 4211, 4280, 4401, 4584))
  expect_equal(sort(out$T_i),
               c(70, 70, 80, 80, 80, 80, 80, 80, 100, 100, 100, 100, 100, 100,
                 100, 100, 100, 100, 120, 120, 120, 120, 120, 120, 120, 120, 120,
                 120, 120, 120))
  expect_equal(sort(out$sim_i),
               c(0.604150050986692, 0.6052739242771, 0.610176794203946, 0.612023739225788,
                 0.614492285215261, 0.620227553733224, 0.621164966971297, 0.622783036449866,
                 0.627319873312495, 0.633248513173231, 0.635107424444073, 0.635521689774035,
                 0.636802853348891, 0.647105239817411, 0.653594865973984, 0.65658700595948,
                 0.661407919988764, 0.66171351823905, 0.669422058643667, 0.677205987228789,
                 0.677569883368504, 0.681830754313774, 0.682351757248779, 0.682661399905778,
                 0.685176833969805, 0.698700404217461, 0.70183423417834, 0.743425396253583,
                 0.755487202273469, 0.764657255410779))

})



################################################################################
################################################################################
################################################################################

context("Testing finetune_maxima")

test_that("Suite 1. Finetune_maxima unchanged", {

  set.seed(20181018)

  x0 <- cos(seq(0, 2 * pi * 5, length.out = 1001))
  finetune.maxima.x  <- x0 + rnorm(1001, sd = 0.3)
  tau.TMP <- 200
  s.TMP <- 200

  nbh.wing <- 30
  x.Fitted <- rep(NA, length(finetune.maxima.x))
  template.vl.min <- 100
  template.vl.max <- 300

  out <- finetune_maxima(s.TMP,
                         tau.TMP,
                         nbh.wing,
                         x.Fitted,
                         finetune.maxima.x,
                         template.vl.min,
                         template.vl.max)

  expect_equal(out, c(212, 193))

})

