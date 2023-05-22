


#' Pattern Segmentation From a Time-series via ADEPT
#'
#' Segment pattern from a time-series \code{x} via Adaptive Empirical Pattern
#' Transformation (ADEPT).
#'
#' @param x A numeric vector. A time-series to segment pattern from.
#' @param x.fs A numeric scalar. Frequency at which a time-series \code{x} is collected,
#' expressed in a number of observations per second.
#' @param template A list of numeric vectors, or a numeric vector.
#' Each vector represents a distinct pattern template used in segmentation.
#' @param pattern.dur.seq  A numeric vector. A grid of potential pattern durations
#' used in segmentation. Expressed in seconds. See: Details.
#' @param similarity.measure A character scalar. Statistic used to compute similarity
#' between a time-series \code{x}  and pattern templates. Currently supported values:
#' \itemize{
#'   \item \code{"cov"} - covariance,
#'   \item \code{"cor"} - correlation,
#' }
#' Default is \code{"cov"}.
#' @param similarity.measure.thresh A numeric scalar. Threshold of minimal similarity
#' value between a time-series \code{x} and a template
#' below which the algorithm does not identify a pattern occurrence from \code{x}.
#' Default is \code{0}.
#' @param x.adept.ma.W A numeric scalar.
#' A length of a window used in moving average smoothing of a time-series \code{x} for
#'   similarity matrix computation. Expressed in seconds.
#'  Default is \code{NULL} (no smoothing applied).
#' @param finetune A character scalar. A type of fine-tuning procedure employed in
#'  segmentation. Defaults to \code{NULL} (no fine-tuning procedure employed). Currently supported values:
#' \itemize{
#'   \item \code{"maxima"} - tunes preliminarily identified
#'   beginning and end of a pattern so
#'   as they correspond to local maxima of time-series \code{x} (or smoothed version of \code{x})
#'   found within neighbourhoods of preliminary locations.
#' }
#' @param finetune.maxima.ma.W A numeric scalar.
#' A length of a window used in moving average smoothing of a time-series \code{x} in
#'  \code{"maxima"} fine-tuning procedure. Expressed in seconds.
#'  Default is \code{NULL} (no smoothing applied).
#' @param finetune.maxima.nbh.W A numeric scalar.
#' A length of the two neighborhoods centered at preliminarily identified
#' beginning and end of a pattern
#' within which we search for local maxima of \code{x} (or smoothed version of \code{x}) in \code{"maxima"}
#' fine-tuning procedure. Expressed in seconds. Default is \code{NULL}.
#' Note: if the length provided corresponds to an even number of \code{x} vector indices,
#' it will be rounded down so as the corresponding number of vector indices is its closest odd number.
#' @param run.parallel A logical scalar. Whether or not to use parallel execution in the algorithm
#' with \code{parallel} package. Default is \code{FALSE}. DOES NOT WORK ON WINDOWS.
#' @param run.parallel.cores An integer scalar.
#' The number of cores to use for parallel execution. Defaults to 1L (no parallel).
#' DOES NOT WORK ON WINDOWS.
#' @param x.cut  A logical scalar. Whether or not to use time optimization procedure in
#' which a time-series \code{x} is cut into parts and segmentation is performed for
#' each part of \code{x} separately. Recommended for a time-series \code{x} of vector length
#'  above 30,000. Default is \code{TRUE}.
#' @param x.cut.vl An integer scalar.
#' Defines a vector length of parts that \code{x} vector is cut into during the execution time optimization procedure.
#' Default is \code{6000} (recommended).
#'
#' @details
#' Function implements Adaptive Empirical Pattern Transformation (ADEPT) method for pattern segmentation
#' from a time-series \code{x}.
#' ADEPT is optimized to perform fast, accurate walking strides segmentation from
#' high-density data collected with a wearable accelerometer during walking.
#'
#'  ADEPT identifies patterns in a time-series \code{x} via maximization of chosen
#'  similarity statistic (correlation, covariance, etc.) between a time-series \code{x}
#'  and a pattern template(s). It accounts for variability in both
#'  (1) pattern duration and (2) pattern shape.
#'
#' @return A \code{data.frame} with segmentation results. Each row
#' describes one identified pattern occurrence:
#' \itemize{
#'   \item \code{tau_i} - index of  \code{x} where pattern starts,
#'   \item \code{T_i} - pattern duration, expressed in \code{x} vector length,
#'   \item \code{sim_i} -  similarity between a pattern and \code{x};
#'   note: if \code{"maxima"} fine-tune and/or \code{x} smoothing is employed,
#'   the similarity value between the final segmented pattern and a template
#'    may differ from the value in this table,
#'   \item \code{template_i} - index of a template best matched to \code{x}.
#' }
#'
#' @export
#'
#' @importFrom dplyr arrange mutate lag filter select
#' @importFrom magrittr '%>%'
#' @importFrom parallel mclapply
#'
#' @references
#' Karas, M., Straczkiewicz, M., Fadel, W., Harezlak, J., Crainiceanu, C.M.,
#' Urbanek, J.K. (2019). Adaptive empirical pattern
#' transformation (ADEPT) with application to walking stride segmentation.
#' Biostatistics. https://doi.org/10.1093/biostatistics/kxz033
#'
#' @examples
#' ## Example 1: Simulate a time-series `x`. Assume that
#' ## - `x` is collected at a frequency of 100 Hz,
#' ## - there is one shape of pattern present within `x`,
#' ## - each pattern lasts 1 second,
#' ## - there is no noise in the collected data.
#' true.pattern <- cos(seq(0, 2 * pi, length.out = 100))
#' x <- c(true.pattern[1], replicate(10, true.pattern[-1]))
#' ## Segment pattern from x.
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = true.pattern,
#'   pattern.dur.seq = c(0.9, 0.95, 1.03, 1.1),
#'   similarity.measure = "cor")
#' out
#' ## Segment pattern from x. Now assume a grid of potential pattern duratios
#' ## contains true pattern duration
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = true.pattern,
#'   pattern.dur.seq = c(0.9, 0.95, 1, 1.03, 1.1),
#'   similarity.measure = "cor")
#' out
#'
#' ## Example 2: Simulate a time-series `x`. Assume that
#' ## - `x` is collected at a frequency of 100 Hz,
#' ## - there are two shapes of pattern present within `x`,
#' ## - patterns have various duration,
#' ## - there is no noise in the collected data.
#' true.pattern.1 <- cos(seq(0, 2 * pi, length.out = 200))
#' true.pattern.2 <- true.pattern.1
#' true.pattern.2[70:130] <- 2 * true.pattern.2[min(70:130)] + abs(true.pattern.2[70:130])
#' x <- numeric()
#' for (vl in seq(70, 130, by = 10)){
#'   true.pattern.1.s <- approx(
#'     seq(0, 1, length.out = 200),
#'     true.pattern.1, xout = seq(0, 1, length.out = vl))$y
#'   true.pattern.2.s <- approx(
#'     seq(0, 1, length.out = 200),
#'     true.pattern.2, xout = seq(0, 1, length.out = vl))$y
#'   x <- c(x, true.pattern.1.s[-1], true.pattern.2.s[-1])
#'   if (vl == 70) x <- c(true.pattern.1.s[1], x)
#' }
#' ## Segment pattern from x. Use a `template` object consisting of both
#' ## true patterns used in `x` simulation.
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = list(true.pattern.1, true.pattern.2),
#'   pattern.dur.seq = 60:130 * 0.01,
#'   similarity.measure = "cor")
#' out
#'
#' ## Example 3: Simulate a time-series `x`. Assume that
#' ## - `x` is collected at a frequency of 100 Hz,
#' ## - there are two shapes of a pattern present within `x`,
#' ## - patterns have various duration,
#' ## - there is noise in the collected data.
#' set.seed(1)
#' x <- x + rnorm(length(x), sd = 0.5)
#' ## Segment pattern from x.
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = list(true.pattern.1, true.pattern.2),
#'   pattern.dur.seq =  60:130 * 0.01,
#'   similarity.measure = "cor")
#' out
#' ## Segment pattern from x. Use `x.adept.ma.W` to define a length of a smoothing
#' ## window to smooth `x` for similarity matrix computation.
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = list(true.pattern.1, true.pattern.2),
#'   pattern.dur.seq =  60:130 * 0.01,
#'   similarity.measure = "cor",
#'   x.adept.ma.W = 0.1)
#' out
#' ## Segment pattern from x. Use `x.adept.ma.W` to define a length of a smoothing
#' ## window to smooth `x` for similarity matrix computation. Employ a fine-tuning
#' ## procedure for stride identification.
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = list(true.pattern.1, true.pattern.2),
#'   pattern.dur.seq =  60:130 * 0.01,
#'   similarity.measure = "cor",
#'   x.adept.ma.W = 0.1,
#'   finetune = "maxima",
#'   finetune.maxima.nbh.W = 0.3)
#' out
#' ## Segment pattern from x. Employ a fine-tuning procedure for stride
#' ## identification. Smooth `x` for both similarity matrix computation
#' ## (set `x.adept.ma.W = 0.1`) and for  fine-tune peak detection procedure
#' ## (set `finetune.maxima.nbh.W = 0.3`).
#' out <- segmentPattern(
#'   x = x,
#'   x.fs = 100,
#'   template = list(true.pattern.1, true.pattern.2),
#'   pattern.dur.seq =  60:130 * 0.01,
#'   similarity.measure = "cor",
#'   x.adept.ma.W = 0.1,
#'   finetune = "maxima",
#'   finetune.maxima.nbh.W = 0.3)
#' out
#'
segmentPattern <- function(x,
                           x.fs,
                           template,
                           pattern.dur.seq,
                           similarity.measure = "cov",
                           similarity.measure.thresh = 0.0,
                           x.adept.ma.W = NULL,
                           finetune = NULL,
                           finetune.maxima.ma.W = NULL,
                           finetune.maxima.nbh.W = NULL,
                           run.parallel = FALSE,
                           run.parallel.cores = 1L,
                           x.cut = TRUE,
                           x.cut.vl = 6000){


  ## ---------------------------------------------------------------------------
  ## Check if correct objects were passed to the function
  x.cut.vl <- as.integer(x.cut.vl)
  if(!is.null(run.parallel.cores)) run.parallel.cores <-  as.integer(run.parallel.cores)
  if (!(all(is.numeric(x)) & is.atomic(x))) stop("x must be a numeric (atomic) vector.")
  if (!(length(x.fs) == 1 & is.numeric(x.fs) & x.fs > 0 & is.atomic(x.fs))) stop("x.fs must be a positive numeric scalar.")
  template.cond1 <- all(is.numeric(template)) & is.atomic(template)
  template.cond2 <- is.list(template) & all(sapply(template, function(vec) all(is.numeric(vec)) & is.atomic(vec)))
  if (!(template.cond1 || template.cond2)) stop("template must be a numeric (atomic) vector, or a list of numeric (atomic) vectors.")
  if (!(all(is.numeric(pattern.dur.seq)) & is.atomic(pattern.dur.seq) & all(pattern.dur.seq > 0))) stop("pattern.dur.seq must be a numeric (atomic) vector of positive values.")
  if (!(similarity.measure %in% c("cov", "cor"))) stop("similarity.measure must be one of: 'cov', 'cor'.")
  if (!(is.null(x.adept.ma.W) || (length(x.adept.ma.W) == 1 & is.numeric(x.adept.ma.W) & x.adept.ma.W > 0))) stop("x.adept.ma.W must be NULL or a positive numeric scalar.")
  if (!(is.null(finetune) || finetune == "maxima")) stop("finetune must be NULL or 'maxima'.")
  if (!(is.null(finetune.maxima.ma.W) || (length(finetune.maxima.ma.W) == 1 & is.numeric(finetune.maxima.ma.W) & finetune.maxima.ma.W > 0))) stop("finetune.maxima.ma.W must be NULL or a positive numeric scalar.")
  if (!(is.null(finetune.maxima.nbh.W) || (length(finetune.maxima.nbh.W) == 1 & is.numeric(finetune.maxima.nbh.W) & finetune.maxima.nbh.W > 0))) stop("finetune.maxima.nbh.W must be NULL or a positive numeric scalar.")
  if (!(run.parallel %in% c(TRUE, FALSE))) stop("run.parallel must be a logical scalar.")
  if (!(is.null(run.parallel.cores) || (length(run.parallel.cores) == 1 & is.integer(run.parallel.cores) & run.parallel.cores > 0))) stop("run.parallel.cores must me NULL or a positive integer scalar")
  if (!(length(x.cut) == 1 & x.cut %in% c(TRUE, FALSE))) stop("x.cut must be a logical scalar.")
  if (!(is.null(x.cut.vl) || (length(x.cut.vl) == 1 & is.integer(x.cut.vl) & x.cut.vl > 0))) stop("x.cut.vl must me NULL or a positive integer scalar")


  ## ---------------------------------------------------------------------------
  ## Compute a list of rescaled template(s)

  ## Define grid of template vector lengths (corresponding to pattern durations)
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
  template.vl.max <- max(template.vl)
  template.vl.min <- min(template.vl)
  ## Rescale templates
  if (!is.list(template)) template <- list(template)
  template.scaled <- scaleTemplate(template, template.vl)


  ## ---------------------------------------------------------------------------
  ## Smooth x signal for ADEPT similarity matrix computation

  if (!is.null(x.adept.ma.W)){
    # W.vl       <- x.adept.ma.W * x.fs
    x.smoothed <- get.x.smoothed(x = x,
                                 W = x.adept.ma.W,
                                 x.fs = x.fs)
  } else {
    x.smoothed <- x
  }

  ## ---------------------------------------------------------------------------
  ## Smooth x signal for fine tuning

  ## Fine-tuning type-specific procedures
  if (!is.null(finetune) && finetune == "maxima"){

    ## Signal smoothing for fine tunning
    if (!is.null(finetune.maxima.ma.W) && finetune.maxima.ma.W > 0){
      # W.vl              <- finetune.maxima.ma.W * x.fs
      finetune.maxima.x <- get.x.smoothed(x = x,
                                          W = finetune.maxima.ma.W,
                                          x.fs = x.fs)
    } else {
      finetune.maxima.x <- x
    }

    ## Other fine-tuning components
    if (!(finetune.maxima.nbh.W > 0)) stop("finetune.maxima.nbh.W should be greater than 0 for finetune == 'maxima'")
    ## Added round() @MK 2020-01-06
    finetune.maxima.nbh.vl <- round(finetune.maxima.nbh.W * x.fs)
    if (!(finetune.maxima.nbh.vl > 0)) stop("finetune.maxima.nbh.W should be greater")
  }


  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## PER-SEGMENT COMPUTATION

  ## If no signal cutting to parts is allowed
  if (!x.cut) x.cut.vl <- length(x)

  x.cut.margin <- template.vl.max - 1
  # x.cut.margin <- 2 * template.vl.max - 1
  x.cut.seq <- seq(1, to = length(x), by = x.cut.vl)

  template.idx.mat.i <- NULL

  # If you were to run this code on Windows, mclapply would simply call lapply,
  # so the code works but sees no speed gain.
  # source: http://dept.stat.lsa.umich.edu/~jerrick/courses/stat701/notes/parallel.html

  # define number of cores to use in parallel
  mc.cores.val <- ifelse (run.parallel & (!(is.null(run.parallel.cores))), run.parallel.cores, 1L)

  out.list <- parallel::mclapply(x.cut.seq, function(i){
    ## Define current x part indices
    idx.i <- i : min((i + x.cut.vl + x.cut.margin), length(x))
    ## If we cannot fit the longest pattern, return NULL
    if (length(idx.i) <= max(template.vl)) return(NULL)
    ## Compute similarity matrix
    similarity.i <- similarityMatrix(x = x.smoothed[idx.i],
                                     template.scaled = template.scaled,
                                     similarity.measure = similarity.measure)

    similarity.mat.i <- similarity.i$similarity
    template.idx.mat.i <- similarity.i$idx

    ## Run max and tine procedure
    out.df.i <- maxAndTune(x = x[idx.i],
                           template.vl = template.vl,
                           similarity.mat = similarity.mat.i,
                           similarity.measure.thresh = similarity.measure.thresh,
                           template.idx.mat = template.idx.mat.i,
                           finetune = finetune,
                           finetune.maxima.x = finetune.maxima.x[idx.i],
                           finetune.maxima.nbh.vl = finetune.maxima.nbh.vl)

    ## Shift \tau parameter according to which part of signal x we are currently working with
    if (nrow(out.df.i) > 0){
      out.df.i$tau_i <- out.df.i$tau_i + i - 1
      return(out.df.i)

    } else {
      ## Return empty data frame
      return(data.frame(tau_i = numeric(),
                        T_i = numeric(),
                        sim_i = numeric(),
                        template_i = numeric()))
    }
  }, mc.cores = getOption("mc.cores", mc.cores.val))

  ## ---------------------------------------------------------------------------
  ## Clear up after possibly multiple stride occurrences
  out.df <- do.call("rbind", out.list)

  ## To surpress the "Note" on package check
  ## (after https://github.com/Rdatatable/data.table/issues/850)
  tau_i <- NULL; T_i <- NULL; tau_i_diff <- NULL; sim_i <- NULL

  k <- floor((template.vl.max-1)/template.vl.min)
  # k <- floor((2 * template.vl.max - 1)/template.vl.min)
  if (k > 0){
    for (i in 1:k){
      out.df <-
        out.df %>%
        # arrange(tau_i) %>%
        arrange(tau_i, dplyr::desc(sim_i)) %>%
        mutate(tau_i_diff = lag(tau_i + T_i - 1) - tau_i) %>%
        filter(tau_i_diff <= 0 | is.na(tau_i_diff))
    }
    out.df <-
      out.df %>%
      select(-tau_i_diff)
  }

  return(out.df)
}

