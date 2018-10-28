



#' Title
#'
#' @param x
#' @param x.fs
#' @param template
#' @param pattern.dur.min
#' @param pattern.dur.max
#' @param pattern.dur.step
#' @param pattern.dur.seq
#' @param similarity.measure
#' @param similarity.measure.thresh
#' @param finetune
#' @param finetune.maxima.ma.W
#' @param finetune.maxima.nbh.W
#' @param run.parallel
#' @param run.parallel.ncores
#'
#' @return Data frame with summary of segmented pattern occurences; has the
#' following columns:
#' \describe{
#'   \item{One}{First item}
#'   \item{Two}{Second item}
#' }
#' @export
#'
#' @examples
segmentPattern <- function(x,
                           x.fs,
                           template,
                           # pattern.dur.min = NULL,
                           # pattern.dur.max  = NULL,
                           # pattern.dur.step = NULL,
                           pattern.dur.seq = NULL,
                           similarity.measure = "cov",
                           similarity.measure.thresh = 0.0,
                           x.adept.ma.W = NULL,
                           finetune = NULL,
                           finetune.maxima.ma.W = NULL,
                           finetune.maxima.nbh.W = NULL,
                           run.parallel = FALSE,
                           run.parallel.ncores = NULL){

  ## ---------------------------------------------------------------------------
  ## Compute collection of rescaled template(s)

  # ## Define grid of pattern durations, expressed in time [sec]
  # if (is.null(pattern.dur.seq)){
  #   pattern.dur.seq <- seq(pattern.dur.min, pattern.dur.max, by = pattern.dur.step)
  # }

  ## Define grid of template vector lengths (corresponding to pattern durations)
  template.vl <- pattern.dur.seq * x.fs
  template.vl <- sort(unique(round(template.vl)))
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



  ###############################################
  ###   THIS BELOW MAY be x SUBSET-SPECIFIC   ###
  ###############################################


  ## ---------------------------------------------------------------------------
  ## ADEPT: covariance matrix computation

  ## Compute similarity matrix
  similarity.mat <- adeptSimilarity(x.smoothed,
                                    template.scaled,
                                    similarity.measure,
                                    run.parallel,
                                    run.parallel.ncores)

  ## Sanity check
  if (!(ncol(similarity.mat) == length(x))) stop("!(ncol(similarity.mat) == length(x))")
  if (!(nrow(similarity.mat) == length(template.vl))) stop("!(nrow(similarity.mat) == length(template.vl))")

  ## ---------------------------------------------------------------------------
  ## ADEPT: maximization tuning algorithm

  out.df <- maxAndTune(x,
                       template.vl,
                       similarity.mat,
                       similarity.measure.thresh,
                       finetune,
                       finetune.maxima.x,
                       finetune.maxima.nbh.vl)


  ## ---------------------------------------------------------------------------
  ## return the result
  out.df <- out.df[order(out.df$tau_i), ]
  return(out.df)

}




