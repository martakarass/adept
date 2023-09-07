run_both = function(func = segmentPattern, args) {
  args$fast = FALSE
  time_out = system.time({
    out = do.call(func, args = args)
  })
  args$fast = TRUE
  time_out_fast = system.time({
    out_fast = do.call(func, args = args)
  })
  list(
    out = out,
    out_fast = out_fast,
    time_out = time_out,
    time_out_fast = time_out_fast
  )
}
test_both = function(func = segmentPattern, args) {
  res = run_both(func = func, args = args)
  testthat::expect_equal(res$out, res$out_fast)
}

testthat::test_that("fast segmentPattern gives same results", {
  ## Example 1: Simulate a time-series `x`. Assume that
  ## - `x` is collected at a frequency of 100 Hz,
  ## - there is one shape of pattern present within `x`,
  ## - each pattern lasts 1 second,
  ## - there is no noise in the collected data.
  true.pattern <- cos(seq(0, 2 * pi, length.out = 100))
  x <- c(true.pattern[1], replicate(10, true.pattern[-1]))
  ## Segment pattern from x.
  args = list(
    x = x,
    x.fs = 100,
    template = true.pattern,
    pattern.dur.seq = c(0.9, 0.95, 1.03, 1.1),
    similarity.measure = "cor"
  )
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

  ## Segment pattern from x. Now assume a grid of potential pattern duratios
  ## contains true pattern duration
  args <- list(
    x = x,
    x.fs = 100,
    template = true.pattern,
    pattern.dur.seq = c(0.9, 0.95, 1, 1.03, 1.1),
    similarity.measure = "cor",
    compute.template.idx = FALSE)
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

  ## Example 2: Simulate a time-series `x`. Assume that
  ## - `x` is collected at a frequency of 100 Hz,
  ## - there are two shapes of pattern present within `x`,
  ## - patterns have various duration,
  ## - there is no noise in the collected data.
  true.pattern.1 <- cos(seq(0, 2 * pi, length.out = 200))
  true.pattern.2 <- true.pattern.1
  true.pattern.2[70:130] <- 2 * true.pattern.2[min(70:130)] + abs(true.pattern.2[70:130])
  x <- numeric()
  for (vl in seq(70, 130, by = 10)){
    true.pattern.1.s <- approx(
      seq(0, 1, length.out = 200),
      true.pattern.1, xout = seq(0, 1, length.out = vl))$y
    true.pattern.2.s <- approx(
      seq(0, 1, length.out = 200),
      true.pattern.2, xout = seq(0, 1, length.out = vl))$y
    x <- c(x, true.pattern.1.s[-1], true.pattern.2.s[-1])
    if (vl == 70) x <- c(true.pattern.1.s[1], x)
  }
  ## Segment pattern from x. Use a `template` object consisting of both
  ## true patterns used in `x` simulation.
  args <- list(
    x = x,
    x.fs = 100,
    template = list(true.pattern.1, true.pattern.2),
    pattern.dur.seq = 60:130 * 0.01,
    similarity.measure = "cor",
    compute.template.idx = FALSE)
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

  ## Example 3: Simulate a time-series `x`. Assume that
  ## - `x` is collected at a frequency of 100 Hz,
  ## - there are two shapes of a pattern present within `x`,
  ## - patterns have various duration,
  ## - there is noise in the collected data.
  set.seed(1)
  x <- x + rnorm(length(x), sd = 0.5)
  ## Segment pattern from x.
  args <- list(
    x = x,
    x.fs = 100,
    template = list(true.pattern.1, true.pattern.2),
    pattern.dur.seq =  60:130 * 0.01,
    similarity.measure = "cor",
    compute.template.idx = FALSE)
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

  ## Segment pattern from x. Use `x.adept.ma.W` to define a length of a smoothing
  ## window to smooth `x` for similarity matrix computation.
  args <- list(
    x = x,
    x.fs = 100,
    template = list(true.pattern.1, true.pattern.2),
    pattern.dur.seq =  60:130 * 0.01,
    similarity.measure = "cor",
    x.adept.ma.W = 0.1,
    compute.template.idx = FALSE)
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

  ## Segment pattern from x. Use `x.adept.ma.W` to define a length of a smoothing
  ## window to smooth `x` for similarity matrix computation. Employ a fine-tuning
  ## procedure for stride identification.
  args <- list(
    x = x,
    x.fs = 100,
    template = list(true.pattern.1, true.pattern.2),
    pattern.dur.seq =  60:130 * 0.01,
    similarity.measure = "cor",
    x.adept.ma.W = 0.1,
    finetune = "maxima",
    finetune.maxima.nbh.W = 0.3,
    compute.template.idx = FALSE)
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

  ## Segment pattern from x. Employ a fine-tuning procedure for stride
  ## identification. Smooth `x` for both similarity matrix computation
  ## (set `x.adept.ma.W = 0.1`) and for  fine-tune peak detection procedure
  ## (set `finetune.maxima.nbh.W = 0.3`).
  args <- list(
    x = x,
    x.fs = 100,
    template = list(true.pattern.1, true.pattern.2),
    pattern.dur.seq =  60:130 * 0.01,
    similarity.measure = "cor",
    x.adept.ma.W = 0.1,
    finetune = "maxima",
    finetune.maxima.nbh.W = 0.3,
    compute.template.idx = FALSE)
  res = test_both(args = args)
  args$similarity.measure = "cov"
  res = test_both(args = args)

})


testthat::test_that("fast segmentWalking gives same result", {
  testthat::skip_if_not_installed("adeptdata")
  library(dplyr)
  xyz <- adeptdata::acc_walking_IU %>%
    filter(loc_id == "left_wrist", subj_id == "id86237981") %>%
    arrange(time_s) %>%
    select(v1 = x, v2 = y, v3 = z) %>%
    as.matrix()
  # define raw accelerometry data sample frequency
  xyz.fs <- 100

  # define template list based on predefined templates
  template_mat <- adeptdata::stride_template$left_wrist[[3]]
  template <- list(
    template_mat[1, ],
    template_mat[2, ],
    template_mat[3, ]
  )

  # run walking segmentation
  # (parallel supported, except for Windows; see run.parallel, run.parallel.cores args)
  args = list(xyz = xyz,
              xyz.fs = xyz.fs,
              template = template)
  test_both(func = segmentWalking, args = args)
  res = run_both(func = segmentWalking, args = args)

})
