


#' Pattern Segmentation From a Time-series via ADEPT
#'
#' Segment pattern from a time-series \code{x} via Adaptive Empirical Pattern
#' Transformation (ADEPT).
#'
#' @param x A numeric vector. A time-series to segment pattern occurrences from.
#' @param x.fs A numeric scalar. Frequency at which a time-series \code{x} is collected,
#' expressed in a number of observations per second.
#' @param template A list of numeric vectors, or a numeric vector.
#' Each vector represents a distinct pattern template used in segmentation.
#' @param pattern.dur.seq  A numeric vector. A grid of pattern duration
#' times used in segmentation. Expressed in seconds. See: Details.
#' @param similarity.measure A character scalar. Statistic used to compute similarity
#' between a time-series \code{x}  and pattern templates. Currently supported values:
#' \itemize{
#'   \item \code{"cov"} - covariance,
#'   \item \code{"cor"} - correlation,
#' }
#' Default is \code{"cov"}.
#' @param similarity.measure.thresh A numeric scalar. Threshold of minimal similarity
#' value between a time-series \code{x} and pattern templates
#' below which the algorithm does not identify a pattern occurrence.
#' Default is \code{0}.
#' @param x.adept.ma.W A numeric scalar.
#' A length of a window used in moving average smoothing of a time-series \code{x} for
#'   similarity matrix computation. Expressed in seconds.
#'  Default is \code{NULL} (no smoothing applied).
#' @param finetune A character scalar. A type of fine-tuning procedure employed in
#'  segmentation. Defaults to \code{NULL} (no fine-tuning procedure employed). Currently supported values:
#' \itemize{
#'   \item \code{"maxima"} - tunes preliminarily identified locations of pattern occurence
#'   beginning and end so
#'   as they correspond to local maxima of time-series \code{x} (or smoothed version of \code{x})
#'   found within neighbourhoods of preliminary locations.
#' }
#' @param finetune.maxima.ma.W A numeric scalar.
#' A length of a window used in moving average smoothing of a time-series \code{x} in
#'  \code{"maxima"} fine-tuning procedure. Expressed in seconds.
#'  Default is \code{NULL} (no smoothing applied).
#' @param finetune.maxima.nbh.W A numeric scalar.
#' A length of the two neighborhoods centered at preliminarily identified pattern occurrence beginning and end points
#' within which we search for local maxima of \code{x} (or smoothed version of \code{x}) in \code{"maxima"}
#' fine-tuning procedure. Expressed in seconds. Default is \code{NULL}.
#' Note: if the length provided corresponds to an even number of \code{x} vector indices,
#' it will be rounded down so as the corresponding number of vector indices is its closest odd number.
#' @param run.parallel A logical scalar. Whether or not to use parallel execution in the algorithm.
#' The \code{future} package
#' is used to execute code asynchronously. Default is \code{FALSE}.
#' @param run.parallel.cores An integer scalar.
#' The number of cores to use for parallel execution.
#' Default is \code{NULL}. If not specified, the number of cores is set to a number of
#' cores available minus 1.
#' @param x.cut  A logical scalar. Whether or not to use time optimization procedure in
#' which a time-series \code{x} is cut into parts and segmentation is performed for
#' each part of \code{x} separately. Recommended for a time-series \code{x} of vector length
#'  above 30,000. Default is \code{TRUE}.
#' @param x.cut.vl An integer scalar.
#' Defines a vector length of parts that \code{x} vector is cut into during the execution time optimization procedure.
#' Default is \code{6000} (recommended).
#' @param compute.template.idx A logical scalar. Whether or not to compute and return information about
#' which of the provided pattern templates yielded a similarity matrix value
#' that corresponds to an identified pattern occurrence.
#' Setting to \code{TRUE} may increase computation time. Default is \code{FALSE}.
#'
#' @details
#'     Function implements Adaptive Empirical Pattern Transformation (ADEPT) method for pattern segmentation
#' from a time-series \code{x}.
#' ADEPT was designed with the aim of performing fast, accurate walking strides segmentation
#' from high-density data
#' collected from wearable accelerometer worn during continuous walking activity.
#'
#'     ADEPT identifies pattern occurrences from a time-series \code{x} via maximizing similarity
#' (correlation, covariance etc.) between a time-series \code{x} and pattern
#' templates.
#' \itemize{
#'   \item Pattern template is scaled to various scale parameters to allow
#'   for a potentially better match with pattern occurrences that may change
#'   its duration over time.  \code{pattern.dur.seq} argument defines a grid of a
#'    pattern duration to be used in segmentation. Note:  the more dense grid may potentially
#'    increase segmentation accuracy but may also increase method execution time.
#'   \item Multiple pattern templates are allowed simultaneously to
#'   allow for potentially better match with pattern occurrences
#'   that may change its shape over time.
#' }
#' In practice, a pre-defined pattern template may be derived as an empirical pattern, that is,
#' a data-derived vector representing a pattern of interest.
#'
#' @return A \code{data.frame} with segmentation results. Each row
#' of the returned \code{data.frame} describes one identified pattern occurrence:
#' \itemize{
#'   \item \code{tau_i} - index of a time-series \code{x} where identified pattern occurence starts,
#'   \item \code{T_i} - duration of identified pattern occurence, expressed in  a time-series \code{x} vector length,
#'   \item \code{sim_i} -  value of similarity statistic between an identified pattern occurence and corresponding
#'   window of a time-series used in similarity matrix computation;
#'   note: this value corresponds to similarity statistic between
#'   preliminarily identified pattern occurence and corresponding window of a time-series used in similarity matrix computation;
#'   specifically: if the fine-tune procedure is employed,
#'   the similarity value between the final pattern occurence location and corresponding window of time-series \code{x}
#'   signal may differ from the value in this table,
#'   \item \code{template_i} - if \code{compute.template.idx} equals \code{TRUE}:
#'   index of pattern template that yielded a similarity matrix value
#'  corresponding to an identified pattern occurrence;
#'  if \code{compute.template.idx} equals \code{FALSE}: \code{NA}.
#' }
#'
#' @export
#'
#' @import future
#' @importFrom dplyr arrange mutate lag filter select
#' @importFrom magrittr '%>%'
#'
#' @examples
#' ## Example 1:
#' ## - no noise in time-series x generation,
#' ## - all pattern occurences of the same length (101) n time-series x generation.
#' ## Generate signal and template.
#' x0 <- cos(seq(0, 2 * pi * 10, length.out = 1001))
#' x  <- x0
#' template <- x0[1:101]
#' ## Include true pattern occurrence length 101
#' ## (and some redundat for the sake of example).
#' pattern.dur.seq <- c(90, 100, 101, 102, 110)
#' ## Use segmentPattern function to identify beginnings tau_i and duration T_i
#' ## of pattern occurences within a signal x.
#' out <- segmentPattern(x = x,
#'                       x.fs = 1,
#'                       template = template,
#'                       pattern.dur.seq = pattern.dur.seq,
#'                       similarity.measure = "cor")
#' out
#'
#' ## Example 2:
#' ## - no noise in time-series x generation,
#' ## - use pattern occurences of different length in time-series x generation.
#' ## Generate signal and template.
#' set.seed(1)
#' ## Grid of different true pattern occurence durations.
#' s.grid <- sample(60:120, size = 10)
#' x_block <- cos(seq(0, 2 * pi, length.out = 200))
#' ## Generate signal x that consists of "glued" pattern occurrences of different length
#' x <- numeric()
#' for (ss in s.grid){
#'   x_block_interpolated <- approx(seq(0, 1, length.out = 200),
#'                                  x_block,
#'                                  xout = seq(0, 1, length.out = ss))$y
#'   if (length(x)>0){
#'     x <- x[-length(x)]
#'   }
#'   x <- c(x, x_block_interpolated)
#' }
#' ## Pattern template used in algorithm
#' template <- x_block
#' ## Assume dense grid of pattern occurrence duration.
#' pattern.dur.seq <- 60:120
#' ## Use segmentPattern function to identify beginnings tau_i and duration T_i.
#' out <- segmentPattern(x = x,
#'                       x.fs = 1,
#'                       template = template,
#'                       pattern.dur.seq = pattern.dur.seq,
#'                       similarity.measure = "cor")
#' out
#'
#' ## Example 3(a):
#' ## - add noise in time-series x generation,
#' ## - use pattern occurences of different length in time-series x generation.
#' ## Generate signal and template
#' s.grid <- sample(60:120, size = 10)
#' x_block <- cos(seq(0, 2 * pi, length.out = 200))
#' set.seed(1)
#' x <- numeric()
#' for (ss in s.grid){
#'   x_block_interpolated <- approx(seq(0, 1, length.out = 200),
#'                                  x_block,
#'                                  xout = seq(0, 1, length.out = ss))$y
#'   if (length(x)>0){
#'     x <- x[-length(x)]
#'   }
#'   x <- c(x, x_block_interpolated)
#' }
#' x <- x + rnorm(length(x), sd = 0.3)
#' pattern.dur.seq <- seq(50, 150, by = 5)
#' ## Pattern template used in algorithm
#' template <- x_block
#' ## Use segmentPattern function to identify beginnings tau_i and duration T_i:
#' ## - use fine-tune "maxima" procedure.
#' out <- segmentPattern(x = x,
#'                       x.fs = 1,
#'                       template = template,
#'                       pattern.dur.seq = pattern.dur.seq,
#'                       similarity.measure = "cor",
#'                       finetune = "maxima",
#'                       finetune.maxima.ma.W = 30,
#'                       finetune.maxima.nbh.W = 120)
#' out
#' ## Example 3(b):
#' ## Use segmentPattern function to identify beginnings tau_i and duration T_i:
#' ## - use time-series x smooting for ADEPT similarity matrix computation,
#' ## - use fine-tune "maxima" procedure.
#' out <- segmentPattern(x = x,
#'                       x.fs = 1,
#'                       template = template,
#'                       pattern.dur.seq = pattern.dur.seq,
#'                       similarity.measure = "cor",
#'                       x.adept.ma.W = 30,
#'                       finetune = "maxima",
#'                       finetune.maxima.ma.W = 30,
#'                       finetune.maxima.nbh.W = 120)
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
                           run.parallel.cores = NULL,
                           x.cut = TRUE,
                           x.cut.vl = 6000,
                           compute.template.idx = FALSE){

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
    finetune.maxima.nbh.vl <- finetune.maxima.nbh.W * x.fs
  }


  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## ---------------------------------------------------------------------------
  ## PER-SEGMENT COMPUTATION

  ## If no signal cutting to parts is allowed
  if (!x.cut) x.cut.vl <- length(x)

  x.cut.margin <- template.vl.max - 1
  x.cut.seq <- seq(1, to = length(x), by = x.cut.vl)

  template.idx.mat.i <- NULL

  if (run.parallel){
    ## multiproces := multicore, if supported, otherwise multisession
    if (is.null(run.parallel.cores)) run.parallel.cores <- availableCores() - 1
    plan(multiprocess, workers = run.parallel.cores)
  } else {
    plan(sequential)
  }


  out.list.f <- lapply(x.cut.seq, function(i){
    future({
      ## Define current x part indices
      idx.i <- i : min((i + x.cut.vl + x.cut.margin), length(x))
      ## If we cannot fit the longest pattern, return NULL
      if (length(idx.i) <= max(template.vl)) return(NULL)
      ## Compute similarity matrix
      similarity.mat.i <- similarityMatrix(x = x.smoothed[idx.i],
                                           template.scaled = template.scaled,
                                           similarity.measure = similarity.measure)
      ## Compute template index matrix
      if (compute.template.idx){
        template.idx.mat.i <- templateIdxMatrix(x = x.smoothed[idx.i],
                                                template.scaled = template.scaled,
                                                similarity.measure = similarity.measure)
      }
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

    })
  })
  out.list <- lapply(out.list.f, value)


  ## ---------------------------------------------------------------------------
  ## Clear up after possibly multiple stride occurences
  out.df <- do.call("rbind", out.list)

  ## To surpress the "Note" on package check
  ## (after https://github.com/Rdatatable/data.table/issues/850)
  tau_i <- NULL; T_i <- NULL; tau_i_diff <- NULL; sim_i <- NULL

  k <- floor((template.vl.max-1)/template.vl.min)
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

