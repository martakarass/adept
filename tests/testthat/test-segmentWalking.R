
context("Testing segmentWalking function.")

# define raw accelerometry data sample with walking activity
xyz <-
  adeptdata::acc_walking_IU %>%
  filter(loc_id == "left_wrist", subj_id == "id86237981", time_s < 70) %>%
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

test_that("Test that the result remains unchanged.", {

  # run walking segmentation
  out <- segmentWalking(xyz, xyz.fs, template)

  # tau_i
  out_tmp <- as.numeric(table(diff(out$tau_i)))
  out_tmp_exp <- c(1, 1, 1, 8, 11, 10, 11, 15, 5, 3, 1)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(table(diff(out$T_i)))
  out_tmp_exp <- c(1, 1, 3, 4, 4, 12, 4, 8, 6, 10, 7, 4, 1, 1, 1)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(summary(out$sim_i))
  out_tmp_exp <- c(0.609683895990804, 0.937583597090984, 0.95016460341074, 0.936577027289593,
                   0.959373815418658, 0.98301945400386)
  expect_equal(out_tmp, out_tmp_exp)
})


test_that("Test that the result remains unchanged when using parallel.", {

  skip_on_os("windows")

  # run walking segmentation
  out <- segmentWalking(xyz, xyz.fs, template, run.parallel = TRUE, run.parallel.cores = 2)

  # tau_i
  out_tmp <- as.numeric(table(diff(out$tau_i)))
  out_tmp_exp <- c(1, 1, 1, 8, 11, 10, 11, 15, 5, 3, 1)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(table(diff(out$T_i)))
  out_tmp_exp <- c(1, 1, 3, 4, 4, 12, 4, 8, 6, 10, 7, 4, 1, 1, 1)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(summary(out$sim_i))
  out_tmp_exp <- c(0.609683895990804, 0.937583597090984, 0.95016460341074, 0.936577027289593,
                   0.959373815418658, 0.98301945400386)
  expect_equal(out_tmp, out_tmp_exp)
})






