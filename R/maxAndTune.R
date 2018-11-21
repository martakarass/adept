



#' Fine-tune preliminary locations of parttern occurence start and end
#'
#' Fine-tune preliminary locations of parttern occurence start and end
#' via identifying singal peaks in respective neighbourhoods of
#' occurence start and end locations.
#'
#' @param s.TMP Index of preliminarily localized parttern occurence start within a signal.
#' @param tau.TMP Vector length of preliminarily localized parttern occurence within a signal.
#' @param nbh.wing Vector length of each side of a neighbourhood.
#' @param x.Fitted Vector of same length as a signal length; contains information
#' whether or not a signal has a pattern occurence already identified at certain index.
#' @param finetune.maxima.x Singal (or its smoothed version) within which we are
#' searching for pattern occurrences.
#' @param template.vl.min Minimal assumed vector length of a pattern occurrence we are searching for.
#' @param template.vl.max Maximal assumed vector length of a pattern occurrence we are searching for.
#'
#' @return 2-element vector of tuned: (1) index of localized parttern occurence start within a signal,
#' (2) vector length of localized parttern occurence within a signal.
#'
finetune_maxima <- function(s.TMP,
                            tau.TMP,
                            nbh.wing,
                            x.Fitted,
                            finetune.maxima.x,
                            template.vl.min,
                            template.vl.max){

  # message("finetune_maxima2")

  ## Define tau1: pattern occurence start index as identified so far
  ## Define tau2: pattern occurence end index as identified so far
  tau1.TMP <- tau.TMP
  tau2.TMP <- tau.TMP + s.TMP - 1
  x.Fitted.vl <- length(x.Fitted)

  ## Define tau1 neighbourhood
  tau1.nbh <- (max(tau1.TMP - nbh.wing, 1)) : (min(tau1.TMP + nbh.wing, tau2.TMP))
  tau1.nbh.trail <- length(tau1.nbh) - which.max(rev(x.Fitted[tau1.nbh] == 1))
  if (length(tau1.nbh.trail) > 0) {
    tau1.nbh <- tau1.nbh[(tau1.nbh.trail + 2) : length(tau1.nbh)]
  }

  ## Define tau2 neighbourhood
  tau2.nbh <- (max(tau2.TMP - nbh.wing, tau1.TMP)) : (min(tau2.TMP + nbh.wing, x.Fitted.vl))
  tau2.nbh.trail <-  which.max(x.Fitted[tau2.nbh] == 1)
  if (length(tau2.nbh.trail) > 0) {
    tau2.nbh <- tau2.nbh[1 : (tau2.nbh.trail - 1)]
  }

  ## Compute matrix of distances between tau2 and tau1 indices
  ## and define if these are egligible given assumed template vector length range
  tau12.mat       <- outer(tau2.nbh, tau1.nbh, FUN = "-") + 1
  tau12.mat.VALID <- (1 * (tau12.mat <= template.vl.max) + 1 * (tau12.mat >= template.vl.min) - 1)

  ## Identify maxima of `finetune.maxima.x` in the two neighbourhods with within egligible indices
  tau1.nbh.x  <- finetune.maxima.x[tau1.nbh]
  tau2.nbh.x  <- finetune.maxima.x[tau2.nbh]
  x.mat       <- outer(tau2.nbh.x, tau1.nbh.x, FUN = "+")
  x.mat.VALID <- x.mat * tau12.mat.VALID
  which.out   <- which(x.mat.VALID == max(x.mat.VALID), arr.ind = TRUE)[1,]

  tau.NEW    <- tau1.nbh[which.out[2]]
  s.NEW      <- tau2.nbh[which.out[1]] - tau.NEW + 1

  return(c(tau.NEW, s.NEW))
}






## THIS IS (x, similarity.mat) SUBINDICES SPECIFIC
#' Title
#'
#' @param x
#' @param template.vl
#' @param similarity.mat
#' @param similarity.measure.thresh
#' @param finetune
#' @param finetune.maxima.x
#' @param finetune.maxima.nbh.vl
#'
#' @return
#' @export
#'
#' @examples
maxAndTune <- function(x,
                       template.vl,
                       similarity.mat,
                       similarity.measure.thresh,
                       finetune,
                       finetune.maxima.x,
                       finetune.maxima.nbh.vl){

  ## Params
  mat.nrow <- nrow(similarity.mat)
  x.vl     <- length(x)
  template.vl.min <- min(template.vl)
  template.vl.max <- max(template.vl)

  ## Data objects to store iteration results
  out.list <- list()
  x.Fitted <- rep(NA, x.vl)

  ## -------------------------------------------------------------------------
  ## Fine-tuning components

  if (!is.null(finetune) && finetune == "maxima"){
    nbh.wing <- floor((finetune.maxima.nbh.vl + (finetune.maxima.nbh.vl %% 2) - 1)/2)
  }
  ## -------------------------------------------------------------------------

  while (TRUE) {

    ## Check condition if any other pattern is left to be discovered
    rl <- rle(is.na(x.Fitted))
    max.empty <- max(rl$lengths[rl$values])
    if (max.empty < template.vl.min){
      break
    }

    ## Determine current maximum value in similarity matrix
    mat.MAX <- max(similarity.mat, na.rm = TRUE)
    if (mat.MAX < similarity.measure.thresh) {
      break
    }

    ## Identify parameters s and tau corresponding to maximum of covariance matrix
    ## s:   expressed as vector length
    ## tau: expressed as index of x vector
    mat.MAX.IDX <- which(similarity.mat == mat.MAX, arr.ind = TRUE)[1, ]
    tau.TMP     <- mat.MAX.IDX[2]
    s.TMP       <- template.vl[mat.MAX.IDX[1]]


    ## -------------------------------------------------------------------------
    ## Fine-tuning

    if (!is.null(finetune) && finetune == "maxima"){
      finetune.out <- finetune_maxima(s.TMP, tau.TMP, nbh.wing, x.Fitted, finetune.maxima.x,
                                       template.vl.min, template.vl.max)
      tau.TMP <- finetune.out[1]
      s.TMP   <- finetune.out[2]
    }
    ## -------------------------------------------------------------------------


    ## Fill similarity matrix with NA's at locations populated by an identified pattern
    for (i in 1:mat.nrow){
      s.i <-  template.vl[i]
      NArepl.cols.min <- tau.TMP - s.i + 2
      NArepl.cols.min <- min(max(1, NArepl.cols.min), x.vl)
      NArepl.cols.max <- tau.TMP + s.TMP - 2
      NArepl.cols.max <- min(max(1, NArepl.cols.max), x.vl)
      NArepl.cols     <- NArepl.cols.min:NArepl.cols.max
      # print(NArepl.cols)
      similarity.mat[i, NArepl.cols] <- NA
    }
    ## Fill x.Fitted with TRUE at locations populated by an identified pattern
    ## (note we intentionally fit WITHOUT the "joining" points, that is, border points)
    x.Fitted[(tau.TMP + 1):(tau.TMP + s.TMP - 2)] <- 1

    ## Store current iteration-specific results
    out.list[[length(out.list) + 1]] <- c(tau.TMP, s.TMP, mat.MAX)

  }

  ## List to data frame
  out.df <- as.data.frame(do.call(rbind, out.list))
  if (nrow(out.df) > 0) names(out.df) <- c("tau_i", "T_i", "sim_i")

  return(out.df)

}

