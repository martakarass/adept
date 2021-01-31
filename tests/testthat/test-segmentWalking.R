
context("Testing segmentWalking function.")

# define raw accelerometry data sample with walking activity
xyz <-
  adeptdata::acc_walking_IU %>%
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

test_that("Test that the result remains unchanged.", {

  # run walking segmentation
  out <- segmentWalking(xyz, xyz.fs, template)

  # tau_i
  out_tmp <- as.numeric(table(diff(out$tau_i)))
  out_tmp_exp <- c(1, 1, 2, 2, 4, 4, 1, 2, 1, 5, 5, 14, 29, 35, 33, 19, 13, 9,
                 1, 1, 1, 1, 1, 3, 3, 2, 4, 1)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(table(diff(out$T_i)))
  out_tmp_exp <- c(1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1,
                   2, 6, 10, 4, 22, 12, 21, 17, 21, 14, 11, 4, 3, 1, 3, 1, 3, 5,
                   4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(summary(out$sim_i))
  out_tmp_exp <- c(0.447788834120344, 0.927338638590566, 0.950041534695681, 0.936063202907588,
                   0.964910021616661, 0.9871113293843)
  expect_equal(out_tmp, out_tmp_exp)
})


test_that("Test that the result remains unchanged when using parallel.", {

  skip_on_os("windows")

  # run walking segmentation
  out <- segmentWalking(xyz, xyz.fs, template, run.parallel = TRUE, run.parallel.cores = 2)

  # tau_i
  out_tmp <- as.numeric(table(diff(out$tau_i)))
  out_tmp_exp <- c(1, 1, 2, 2, 4, 4, 1, 2, 1, 5, 5, 14, 29, 35, 33, 19, 13, 9,
                   1, 1, 1, 1, 1, 3, 3, 2, 4, 1)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(table(diff(out$T_i)))
  out_tmp_exp <- c(1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1,
                   2, 6, 10, 4, 22, 12, 21, 17, 21, 14, 11, 4, 3, 1, 3, 1, 3, 5,
                   4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2)
  expect_equal(out_tmp, out_tmp_exp)

  # T_i
  out_tmp <- as.numeric(summary(out$sim_i))
  out_tmp_exp <- c(0.447788834120344, 0.927338638590566, 0.950041534695681, 0.936063202907588,
                   0.964910021616661, 0.9871113293843)
  expect_equal(out_tmp, out_tmp_exp)
})






