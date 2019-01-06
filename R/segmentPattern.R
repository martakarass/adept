



#' Pattern Segmentation From a Time-series via ADEPT
#'
#' Perform pattern segmentation from a time-series \code{x} via ADaptive Empirical Pattern
#' Transformation (ADEPT).
#'
#' @details
#' Function implements ADaptive Empirical Pattern Transformation (ADEPT) method for pattern segmentation
#' from a time-series \code{x}.
#' ADEPT was designed with the aim of performing fast, accurate walking strides segmentation from high-density data
#' collected from wearable accelerometer worn during continuous walking activity.
#'
#' ADEPT identifies pattern occurrenes from a time-series \code{x} via maximizing similarity
#' (correlation, covariance etc.) between time-series \code{x} windows and a
#' template vector(s). In practice, a pre-defined pattern template may be derived as an empirical pattern, that is,
#' data-derived vector representing a pattern of interest.
#'
#' To address a possible scenario in which a pattern occurrence is changing its duration
#' within time-series \code{x}, for each window of \code{x} considered, an empirical pattern(s)
#' is scaled to various scale parameters
#' (that is, linearly interpolated into various vector lengths).
#' Also, multiple pattern templates are allowed simultaneously to
#' even better match a pattern which may also change its shape over the time-series \code{x}.
#'
#' @param x A numeric vector; time-series to segment pattern occurrences from.
#' @param x.fs A numeric scalar; frequency at which time-series \code{x} was collected,
#' expressed in a number of observations per second.
#' @param template A numeric vector or a list of numeric vectors. Represents template(s) of a pattern to segment.
#' @param pattern.dur.seq  A numeric vector. Defines a grid of pattern duration time considered in the
#' method. Expressed in seconds. Note: The provided grid is further translated into \code{x} vector length unit;
#' his includes rounding the translated values onto integer values, and (if needed) selecting unique values only.
#' @param similarity.measure A string. Denotes statistic used to define similarity
#' between time-series \code{x} windows and a pattern template. Currently supported values:
#' \describe{
#'   \item{"cov"}{covariance},
#'   \item{"cor"}{correlation}.
#' }
#' Default is \code{"cov"}.
#' @param similarity.measure.thresh A numeric scalar. Defines threshold of minimal similarity
#' value between time-series \code{x} and a pattern template
#' below which we no longer identify a pattern occurrence.
#' @param finetune A string. Defines type of fine-tuning procedure empolyed in
#'  pattern identification. Defaults to \code{NULL}. Currently supported values:
#' \itemize{
#'   \item "maxima" - the procedure tunes preliminarily identified locations of parttern occurence beginning and end so
#'   as they correspond to local maxima of time-series \code{x} (or smoothed version of \code{x})
#'   found within neighbourhoods of those preliminary locations.
#' }
#' @param finetune.maxima.ma.W A numeric scalar.
#' Defines the length of a window used for moving average smoothing of \code{x} signal in \code{"maxima"}
#' fine-tuning procedure. Expressed in seconds. Default is \code{NULL} (no moving average smoothing applied).
#' @param finetune.maxima.nbh.W A numeric scalar.
#' Defines the length of neighborhoods centered at preliminarily identified pattern occurrence beginning and end points
#' within which we search for local maxima of \code{x} (or smoothed version of \code{x}) in \code{"maxima"}
#' fine-tuning procedure. Expressed in seconds. Default is \code{NULL}.
#' Note: if the length provided rounds to an even number of \code{x} vector indices,
#' it will be rounded down so as the corresponding number of vector indices is its closest odd number
#' (and hence neighbourhood parts - to the left and to the right from the  preliminarilyidentified occurence beginning/end
#' point -
#' are of the same length).
#' @param run.parallel Logical. Whether or not to use parallel method execution. The \code{future} package
#' is used to evaluate expressions asynchronously. Default is \code{FALSE}.
#' @param run.parallel.workers An integer scalar.
#' Number of workers (\code{future} package syntax) used in the parallel execution of the method.
#' Default is \code{NULL}. Note: If \code{run.parallel} is \code{TRUE} and
#' \code{run.parallel.workers} is not provided, the number of workers used is set to a number of
#' workers available minus 1.
#' @param x.cut  Logical.
#' Whether or not to use execution time optimization procedure in which time-series \code{x}
#' is cut into parts so as segmentation is performed for each part of \code{x} separately and results aggregated then.
#' Recommended for \code{x}  vector of length above 30,000.
#' (which corresponds to 5 minutes of data collected at frequency 100 Hz).
#' @param x.cut.vl An integer scalar.
#' Defines a vector length of parts onto which \code{x vector is cut in the execution time optimization procedure.
#' Default is 6000 (recommended).
#'
#' @return A \code{data.frame} object with the segmentation result. Each row
#' of the output corresponds to one identified pattern occurrence:
#' \itemize{
#'   \item \code{tau_i} - index of time-series \code{x} where identified pattern occurence starts,
#'   \item \code{T_i} - duration of identified pattern occurence starts, expressed in \code{x} vector length,
#'   \item \code{sim_i} -  value of similarity statistic between identified pattern occurence and corresponding
#'   window of time-series used in similarity matrix computation;
#'   note: this value corresponds to similarity statistic between
#'   preliminarily identified pattern occurence and corresponding window of time-series used in similarity matrix computation;
#'   specifically: if the fine-tune procedure is employed,
#'   the similarity value between the final pattern occurence location and corresponding window of time-series \code{x}
#'   singal may differ from the value in this table.
#' }
#'
#' @export
#'
#' @import future
#'
#' @examples
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
                           run.parallel.workers = NULL,
                           x.cut = TRUE,
                           x.cut.vl = 6000){

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
    x.smoothed <- get.x.smoothed(x, W = x.adept.ma.W, x.fs = x.fs)
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
      finetune.maxima.x <- get.x.smoothed(x, W = finetune.maxima.ma.W, x.fs = x.fs)
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

  x.cut.seq <- seq(1, to = length(x), by = x.cut.vl)
  x.cut.margin <- template.vl.max - 1

  if (run.parallel){
    ## multiproces := multicore, if supported, otherwise multisession
    if (is.null(run.parallel.workers)) run.parallel.workers <- availableCores() - 1
    plan(multiprocess, workers = run.parallel.workers)
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
      similarity.mat.i <- adeptSimilarity(x.smoothed[idx.i],
                                          template.scaled,
                                          similarity.measure)
      ## Run max and tine procedure
      out.df.i <- maxAndTune(x[idx.i],
                             template.vl,
                             similarity.mat.i,
                             similarity.measure.thresh,
                             finetune,
                             finetune.maxima.x[idx.i],
                             finetune.maxima.nbh.vl)
      ## Shift \tau parameter according to which part of signal x we are currently working with
      if (nrow(out.df.i) > 0){
        out.df.i$tau_i <- out.df.i$tau_i + i - 1
        return(out.df.i)
      } else {
        ## Return empty data frame
        return(data.frame(tau_i = numeric(), T_i = numeric(), sim_i = numeric()))
      }
      # if (nrow(out.df.i)>0) out.df.i$tau_i <- out.df.i$tau_i + i - 1
      # return(out.df.i)
    })
  })
  out.list <- lapply(out.list.f, value)


  ## ---------------------------------------------------------------------------
  ## Clear up after possibly multiple stride occurences
  out.df <- do.call("rbind", out.list)

  k <- floor((template.vl.max-1)/template.vl.min)
  if (k > 0){
    for (i in 1:k){
      out.df <-
        out.df %>%
        arrange(tau_i) %>%
        mutate(tau_i_diff = lag(tau_i + T_i - 1) - tau_i) %>%
        filter(tau_i_diff <= 0 | is.na(tau_i_diff))
    }
    out.df <-
      out.df %>%
      select(-tau_i_diff)
  }

  return(out.df)
}




