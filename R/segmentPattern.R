



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
                           pattern.dur.seq,
                           similarity.measure = "cov",
                           similarity.measure.thresh = 0.0,
                           x.adept.ma.W = NULL,
                           finetune = NULL,
                           finetune.maxima.ma.W = NULL,
                           finetune.maxima.nbh.W = NULL,
                           run.parallel = FALSE,
                           run.parallel.ncores = NULL,
                           x.cut = TRUE,
                           x.cut.vl = 6000){

  ## ---------------------------------------------------------------------------
  ## Compute collection of rescaled template(s)

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

    ## -------------------------------------------------------------------------
    ## PARALLEL COMPUTATION

    ## Define number of cores
    if (is.null(run.parallel.ncores)) run.parallel.ncores <- detectCores() - 1
    cl <- makeCluster(run.parallel.ncores)
    ## Export objects to cluster
    clusterExport(cl,
                  c("x.smoothed", "template.scaled", "similarity.measure",
                    "adeptSimilarity", "maxAndTune", "finetune_maxima",
                    "x", "template.vl", "similarity.measure.thresh",
                    "finetune", "finetune.maxima.x", "finetune.maxima.nbh.vl",
                    "x.cut.vl", "x.cut.margin"),
                  envir = environment())

    ## Run ADEPT procedure on each part of x separately
    out.list <- parLapply(cl, x.cut.seq, function(i){
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
      out.df.i$tau_i <- out.df.i$tau_i + i - 1
      return(out.df.i)
    })
    stopCluster(cl)

  } else {

    ## -------------------------------------------------------------------------
    ## NON-PARALLEL COMPUTATION

    ## Run ADEPT procedure on each part of x separately
    out.list <- lapply(x.cut.seq, function(i){
      ## Define current x part indices
      #   ## >>> TESTING x.cut.margin <<<
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
      out.df.i$tau_i <- out.df.i$tau_i + i - 1
      return(out.df.i)
    })

  }

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




