
#' @noRd
#'
repl_NA <- function(vec, undef_val = 0){
  vec[is.na(vec)] <- undef_val
  return(vec)
}

#' @noRd
#'
RunningDiff1 <- function(x){
  out <- c(diff(x), NA)
  return(out)
}

#' @noRd
#'
runningSum <- function(x, W){
  x_ext <- c(x, x[1:(W-1)])
  x_runningSum <- diff(c(0, cumsum(x_ext)), lag = W)
  return(x_runningSum)
}

#' @noRd
#'
RunningMean <- function(x, W, circular = FALSE){
  if (W > length(x)) stop("W must be smaller or equal to the length of x")
  meanx <- runningSum(x, W)/W
  if (!circular){
    meanx[(length(x) - W + 2) : length(x)] <- NA
  }
  return(meanx)
}


#' Walking Stride Pattern Segmentation from Raw Accelerometry Data via ADEPT
#'
#' Segment walking stride pattern from a raw accelerometry data time-series
#' \code{(x,y,z)} via Adaptive Empirical Pattern Transformation (ADEPT).
#' Default algorithm parameters are optimized for a wrist-worn sensor and
#' were evaluated with data collected in the free-living environment.
#'
#' @param xyz A numeric matrix (or data frame) of n × 3 dimension.
#' Three-dimensional
#' raw accelerometry data time-series; accceleration measurements \code{(x,y,z)}
#' collected
#' along three orthogonal axes by the sensor's accelerometer.
#' @param xyz.fs A numeric scalar. Frequency at which a time-series
#' \code{(x,y,z)} is
#' collected, expressed in a number of observations per second.
#' @param template A list of numeric vectors, or a numeric vector. Distinct
#' pattern template(s) of walking stride.
#' @param sim_MIN  numeric scalar. Minimum value of correlation between pattern
#' template(s) and (r_t)_t vector magnitude of accelerometry data.
#' Default used is 0.85.
#' @param dur_MIN A numeric scalar. Minimum value of a stride duration allowed
#' to be identified. Expressed in seconds. Default used is 0.8.
#' @param dur_MAX A numeric scalar. Maximum value of a stride duration allowed
#' to be identified. Expressed in seconds. Default used is 1.4.
#' @param ptp_r_MIN A numeric scalar. Minimum value of "peak to peak"
#'  difference in (r_t)_t vector magnitude data of a stride.
#'  Default used is 0.2.
#' @param ptp_r_MAX A numeric scalar. Maximum value of "peak to peak"
#'  difference in (r_t)_t vector magnitude data of a stride.
#'  Default used is 2.0
#' @param vmc_r_MIN A numeric scalar. Minimum value of VMC in
#' (r_t)_t vector magnitude data
#' of a stride. Default used is 0.05.
#' @param vmc_r_MAX A numeric scalar. Maximum value of VMC in
#' (r_t)_t vector magnitude data
#' of a stride. Default used is 0.5.
#' @param mean_abs_diff_med_p_MAX A numeric scalar. Maximum value of MAD* of
#' Azimuth (az_)_t median for 3 subsequent valid strides.
#' Here, MAD* stands for mean
#' of 2 absolute differences between 3 subsequent values. Default used is 0.5.
#' @param mean_abs_diff_med_t_MAX A numeric scalar. Maximum value of MAD* of
#' Elevation (el_)_t median for 3 subsequent valid strides.
#' Here, MAD* stands for mean
#' of 2 absolute differences between 3 subsequent values. Default used is 0.2.
#' @param mean_abs_diff_dur_MAX  A numeric scalar. Maximum value of MAD* of
#' duration time for 3 subsequent valid strides. Here, MAD* stands for mean
#' of 2 absolute differences between 3 subsequent values. Default used is 0.2.
#' @param run.parallel A logical scalar. Whether or not to use parallel
#' execution in the algorithm
#' with \code{parallel} package. Default is \code{FALSE}.
#' @param run.parallel.cores An integer scalar.
#' The number of cores to use for parallel execution. Defaults to 1L
#' (no parallel).
#' @param compute.template.idx A logical scalar. Whether or not to compute
#' and return information about
#' which of the provided pattern templates yielded a similarity matrix value
#' that corresponds to an identified pattern occurrence.
#' Setting to \code{TRUE} may increase computation time.
#' Default is \code{FALSE}.
#'
#' @return A \code{data.frame} with segmentation results. Each row
#' describes one identified pattern occurrence:
#' \itemize{
#'   \item \code{tau_i} - row index of \code{xyz} where pattern starts,
#'   \item \code{T_i} - pattern duration, expressed in \code{xyz} vector length,
#'   \item \code{sim_i} -  similarity between a pattern and
#'   best-fit template; see \code{segmentPattern} for details,
#'   \item \code{template_i} - if \code{compute.template.idx} equals \code{TRUE}:
#'   index of a template best matched to \code{x};
#'  if \code{compute.template.idx} equals \code{FALSE}: \code{NA},
#'  \item \code{is_walking_i} -- \code{1} if a pattern is identified as walking
#'  stride; \code{0} otherwise.
#' }
#'
#' @export
#'
#' @importFrom pracma cart2sph
#'
#' @examples
#' \donttest{
#'   # define raw accelerometry data sample with walking activity
#'   xyz <-
#'     adeptdata::acc_walking_IU %>%
#'     filter(loc_id == "left_wrist", subj_id == "id86237981") %>%
#'     arrange(time_s) %>%
#'     select(v1 = x, v2 = y, v3 = z) %>%
#'     as.matrix()
#'   # define raw accelerometry data sample frequency
#'   xyz.fs <- 100
#'
#'   # define template list based on predefined templates
#'   template_mat <- adeptdata::stride_template$left_wrist[[3]]
#'   template <- list(
#'     template_mat[1, ],
#'     template_mat[2, ],
#'     template_mat[3, ]
#'   )
#'
#'   # run walking segmentation
#'   segmentWalking(xyz, xyz.fs, template)
#'   segmentWalking(xyz, xyz.fs, template, run.parallel = TRUE, run.parallel.cores = 10)
#' }
#'
segmentWalking <- function(xyz,
                           xyz.fs,
                           template,
                           sim_MIN = 0.85,
                           dur_MIN = 0.8,
                           dur_MAX = 1.4,
                           ptp_r_MIN = 0.2,
                           ptp_r_MAX = 2.0,
                           vmc_r_MIN = 0.05,
                           vmc_r_MAX = 0.5,
                           mean_abs_diff_med_p_MAX = 0.5,
                           mean_abs_diff_med_t_MAX = 0.2,
                           mean_abs_diff_dur_MAX   = 0.2,
                           compute.template.idx = FALSE,
                           run.parallel = FALSE,
                           run.parallel.cores = 1){

  # compute all spherical
  xyz <- as.matrix(xyz)
  xyzptr <- as.data.frame(cbind(xyz, cart2sph(xyz)))
  vm <- xyzptr[, 6]

  # run adept pattern identification
  out <- segmentPattern(
    x = vm,
    x.fs = xyz.fs,
    template = template,
    pattern.dur.seq = seq(0.5, 4, length.out = 30),
    similarity.measure = "cor",
    similarity.measure.thresh = -2,
    x.adept.ma.W = 0.2,
    finetune = "maxima",
    finetune.maxima.ma.W = NULL,
    finetune.maxima.nbh.W = 0.6,
    run.parallel = run.parallel,
    run.parallel.cores = run.parallel.cores,
    x.cut = TRUE,
    x.cut.vl = 6000,
    compute.template.idx = compute.template.idx)

  # generate detailed summary of ADEPT-identified patterns
  out_desc <- matrix(nrow = nrow(out), ncol = 5)
  for (i in 1:nrow(out)){  # i <- 1
    # i-th identified pattern: ADEPT result
    tau_i <- unlist(out[i, "tau_i"])
    T_i   <- unlist(out[i, "T_i"])
    sim_i <- unlist(out[i, "sim_i"])
    idx_i <- tau_i : (tau_i + T_i - 1)
    # pull i-th identified pattern data current
    xyzptr_stride1 <- xyzptr[idx_i, ]
    # summarize i-th identified pattern data current
    out_desc[i, ] <- c(
      median(xyzptr_stride1[,4]), # "med_p"
      median(xyzptr_stride1[,5]), # "med_t"
      diff(range(xyzptr_stride1[, 6])), # "ptp_r"
      mean(abs(xyzptr_stride1[, 6] - mean(xyzptr_stride1[, 6]))), # vmc_r
      T_i / x.fs # dur
    )
  }
  out_desc <- as.data.frame(out_desc)
  names(out_desc) <- c("med_p", "med_t", "ptp_r", "vmc_r", "dur")

  # combine ADEPT result and detailed summary of ADEPT-identified patterns
  ADEPT_F_data <- cbind(out, out_desc)

  # Step 1: filter ADEPT-derived patterns to keep those that:
  # (a) have a correlation with the best-fit template of at least 0.85;
  # (2) have duration time between 0.8 and 1.4 seconds;
  # (3) have (r_t)_t amplitude and mean amplitude deviation (MAD, or VMC) within
  #     [0.2, 2.0]g and [0.05, 2.0]g, respectively
  cond_step1 <-
    (ADEPT_F_data$sim >= sim_MIN) &
    (ADEPT_F_data$dur >= dur_MIN) &  (ADEPT_F_data$dur <= dur_MAX) &
    (ADEPT_F_data$ptp_r >= ptp_r_MIN) &  (ADEPT_F_data$ptp_r <= ptp_r_MAX) &
    (ADEPT_F_data$vmc_r >= vmc_r_MIN) &  (ADEPT_F_data$vmc_r <= vmc_r_MAX)

  # Step 2: filter ADEPT-derived patterns to keep those that
  # occurred in a sequence of at least three that have similar:
  # - duration,
  # - median azimuth angle (az_)_t,
  # - and median elevation angle (el_)_t,
  med_p_diff1 <- abs(repl_NA(RunningDiff1(ADEPT_F_data$med_p), undef_val = Inf))
  med_t_diff1 <- abs(repl_NA(RunningDiff1(ADEPT_F_data$med_t), undef_val = Inf))
  dur_diff1   <- abs(repl_NA(RunningDiff1(ADEPT_F_data$dur), undef_val = Inf))
  # limit on differences between subsequent pattern properties
  cond_stride_mean3      <- (repl_NA(RunningMean(cond_step1, W = 3))) > 0.9
  cond_med_p_diff1_mean2 <- (repl_NA(RunningMean(med_p_diff1, W = 2), mean_abs_diff_med_p_MAX) < mean_abs_diff_med_p_MAX)
  cond_med_t_diff1_mean2 <- (repl_NA(RunningMean(med_t_diff1, W = 2), mean_abs_diff_med_t_MAX) < mean_abs_diff_med_t_MAX)
  cond_dur_diff1_mean2   <- (repl_NA(RunningMean(dur_diff1,   W = 2), mean_abs_diff_dur_MAX) < mean_abs_diff_dur_MAX)
  # determine which patterns pass all the conditions
  cond2_which <- (cond_stride_mean3 & cond_med_p_diff1_mean2 & cond_med_t_diff1_mean2 & cond_dur_diff1_mean2) * 1
  cond2_which <- unlist(lapply(which(cond2_which == 1), function(val) val + c(0,1,2)))
  cond2_which <- unique(cond2_which)
  cond2_which <- cond2_which[cond2_which <= nrow(ADEPT_F_data)]
  cond2_which <- as.integer(cond2_which)
  cond_step2 <- rep(0, nrow(ADEPT_F_data))
  cond_step2[cond2_which] <- 1
  # add information to ADEPT results
  out$is_walking_i <- cond_step2

  return(out)
}

