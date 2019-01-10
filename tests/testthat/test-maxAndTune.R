


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
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cov")

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_equal(sort(out$tau_i),
               c(1, 101, 201, 302, 401, 501, 601, 701, 801, 901))
  expect_equal(sort(out$T_i),
               c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100))

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
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")

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



test_that("Sute 3. Test if maxAndTune works correctly. Case: no noise in signal, no tunning, one of the pattern durations is perfectly.", {

  set.seed(20181019)

  s.grid <- c(80, 100, 120)
  templ <- sin(seq(0, pi, length.out = 200))

  x <- numeric()
  for (ss in sample(s.grid, size = 50, replace = TRUE)){
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
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")

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
               c(1, 100, 179, 278, 377, 456, 575, 694, 773, 852, 951, 1050,
                 1149, 1248, 1347, 1446, 1565, 1664, 1783, 1862, 1981, 2080, 2179,
                 2258, 2337, 2456, 2555, 2654, 2733, 2812, 2931, 3050, 3129, 3208,
                 3287, 3386, 3485, 3604, 3683, 3762, 3841, 3960, 4059, 4158, 4237,
                 4336, 4435, 4514, 4593, 4712))
  expect_equal(sort(out$T_i),
               c(80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80,
                 80, 80, 80, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 120, 120, 120,
                 120, 120, 120, 120, 120, 120, 120, 120, 120))
  expect_equal(sort(out$sim_i),
               c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1, 1))

})





test_that("Sute 4. Test if maxAndTune works correctly. Case: add noise to signal, no tunning", {

  set.seed(20181019)

  s.grid <- c(80, 100, 120)
  templ <- sin(seq(0, pi, length.out = 200))

  x <- numeric()
  for (ss in sample(s.grid, size = 50, replace = TRUE)){
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
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_equal(sort(out$tau_i),
               c(9, 99, 179, 278, 381, 455, 591, 694, 773, 852, 953, 1066, 1148,
                 1250, 1372, 1442, 1566, 1682, 1782, 1861, 1982, 2085, 2185, 2254,
                 2335, 2454, 2553, 2655, 2740, 2809, 2931, 3050, 3131, 3210, 3286,
                 3385, 3486, 3607, 3681, 3760, 3843, 3946, 4065, 4157, 4236, 4335,
                 4436, 4518, 4589, 4712))
  expect_equal(sort(out$T_i),
               c(70, 70, 70, 70, 70, 70, 70, 80, 80, 80, 80, 80, 80, 80, 80,
                 80, 80, 80, 80, 80, 80, 80, 100, 100, 100, 100, 100, 100, 100,
                 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 120, 120,
                 120, 120, 120, 120, 120, 120, 120, 120))
  expect_equal(sort(out$sim_i),
               c(0.340426599740086, 0.415109370300428, 0.419980330235871, 0.456847060035004,
                 0.497335581762142, 0.504037125843788, 0.507362606353015, 0.513700762224663,
                 0.537476478237868, 0.548828861345714, 0.563714458641447, 0.564503073551733,
                 0.565573194390764, 0.567992326644352, 0.569269321870132, 0.571326077006986,
                 0.572212013377846, 0.582703897066748, 0.591679916153738, 0.596907292238624,
                 0.616397992209749, 0.61810377854496, 0.618353415466538, 0.618830385702796,
                 0.619548647998267, 0.62212688187184, 0.628613344387549, 0.640588692090826,
                 0.642285112343817, 0.647445436101857, 0.647834870584385, 0.648556010619438,
                 0.655780168413006, 0.65590008373686, 0.662367552896298, 0.666992021892172,
                 0.669064329835441, 0.677060689537498, 0.682256478373396, 0.689455969518788,
                 0.697764800831258, 0.705855435300397, 0.710353458857682, 0.712414684438509,
                 0.71668401589209, 0.722784398510534, 0.724943316342911, 0.733234924599527,
                 0.73599982508251, 0.742142600411899))

})



test_that("Sute 5. Test if maxAndTune works correctly. Case: add noise to signal, no tunning", {

  set.seed(20181019)

  s.grid <- c(80, 100, 120)
  templ <- sin(seq(0, pi, length.out = 200))

  x <- numeric()
  for (ss in sample(s.grid, size = 50, replace = TRUE)){
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
  similarity.mat <- similarityMatrix(x.smoothed, template.scaled, "cor")

  ## Output
  out <- maxAndTune(x,
                    template.vl,
                    similarity.mat,
                    similarity.measure.thresh,
                    finetune,
                    finetune.maxima.x,
                    finetune.maxima.nbh.vl)

  expect_equal(sort(out$tau_i),
               c(99, 278, 455, 694, 773, 852, 1148, 1250, 1442, 1566, 1782,
                 1861, 1982, 2085, 2254, 2335, 2454, 2553, 2655, 2931, 3050, 3286,
                 3385, 3681, 3843, 3946, 4157, 4335, 4436, 4712))
  expect_equal(sort(out$T_i),
               c(80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 100, 100, 100, 100,
                 100, 100, 100, 100, 100, 100, 100, 100, 100, 120, 120, 120, 120,
                 120, 120, 120))
  expect_equal(sort(out$sim_i),
               c(0.616397992209749, 0.61810377854496, 0.618353415466538, 0.618830385702796,
                 0.619548647998267, 0.62212688187184, 0.628613344387549, 0.640588692090826,
                 0.642285112343817, 0.647445436101857, 0.647834870584385, 0.648556010619438,
                 0.655780168413006, 0.65590008373686, 0.662367552896298, 0.666992021892172,
                 0.669064329835441, 0.677060689537498, 0.682256478373396, 0.689455969518788,
                 0.697764800831258, 0.705855435300397, 0.710353458857682, 0.712414684438509,
                 0.71668401589209, 0.722784398510534, 0.724943316342911, 0.733234924599527,
                 0.73599982508251, 0.742142600411899))

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

