



#' Fine-tune preliminary locations of parttern occurence start and end
#'
#' Fine-tune preliminary locations of parttern occurence start and end
#' via identifying singal peaks in respective neighbourhoods of
#' occurence start and end locations.
#'
#' @details
#' First, respective neighbourhoods of preliminary locations of parttern
#' occurence start and end, \code{tau1.nbh} and \code{tau2.nbh} respectivelty,
#' are defined. To do this:
#' \itemize{
#'   \item we define a vector, centered at preliminary
#'   start and end points, of length (2x \code{nbh.wing} + 1).
#'   \item we assure these neighbourhoods are not outside the scope of, possibly
#'   smoothed, time-series \code{x} (defined here as \code{finetune.maxima.x} vector),
#'   \item we assure these neighbourhoods do not contain areas of time-series \code{x}
#'   where a pattern occurence has been already identified.
#' }
#'
#' Second, we define a matrix \code{tau12.mat} -- a matrix of distances between
#' indices of all pairs of points from \code{tau1.nbh} and \code{tau2.nbh} neighbourhoods.
#' We then identify these pairs of points which are "valid" in a sense the distance
#' between their indices is within ranges defined as minimal and maximal
#' vector length of a pattern occurrence we are searching for.
#' In other words, a stride which has a start and end to be at any "valid"
#' pair of points from \code{tau1.nbh} and \code{tau2.nbh} neighbourhoods, respectively,
#' has vector length within ranges defined as minimal and maximal
#' vector length of a pattern occurrence we are searching for.
#' Matrix \code{tau12.mat.VALID} has entry equal to \code{1} where a corresponding
#' pair of points from \code{tau1.nbh} and \code{tau2.nbh} neighbourhoods is "valid"
#' and has entry equal to \code{0} otherwise.
#'
#' Third, we compute \code{x.mat} -- matrix whose entries are sums of values of,
#' possibly smoothed, time-series \code{x}, (defined here as \code{finetune.maxima.x} vector),
#' computed for each pairs of points from \code{tau1.nbh} and \code{tau2.nbh} neighbourhoods.
#' Clearly, the highest \code{x.mat} matrix entry corresponds to
#' two points for which  \code{finetune.maxima.x} attains its peaks within
#' \code{tau1.nbh} and \code{tau2.nbh} neighbourhoods.
#'
#' Then, we define \code{x.mat.VALID} which is \code{x.mat} subset to a values
#' whose corresponding indices are "valid", as coded within \code{tau12.mat.VALID}
#' matrix. Finally, we identify the pair of points from \code{tau1.nbh} and
#' \code{tau2.nbh} neighbourhoods which correspond to a maximum entry of
#' \code{x.mat.VALID} matrix.
#'
#'
#' @param s.TMP An integer; an index of a preliminarily idenitifed parttern occurence start
#' within a time-series \code{x}.
#' @param tau.TMP An integer; a vector length of a preliminarily idenitifed parttern occurence
#' within a time-series \code{x}.
#' @param nbh.wing An integer; a vector length of each side of a neighbourhood centered
#' at preliminary locations of parttern occurence start and end.
#' @param x.Fitted A numeric vector; has the same length as time-series \code{x} /
#' \code{finetune.maxima.x};  each element of this vector has either value \code{1}
#' if \code{x} has a pattern occurence already
#' identified at a corresponding index, or value \code{NA} otherwise.
#' @param finetune.maxima.x A numeric vector; either time-series \code{x} within which we are
#' searching for pattern occurrences, or its already
#' smoothed version.
#' @param template.vl.min An integer; minimal vector length of a pattern occurrence we are searching for.
#' @param template.vl.max An integer; maximal vector length of a pattern occurrence we are searching for.
#'
#' @return A 2-element numeric vector that contains:
#' \itemize{
#'   \item at its 1st element: an index of localized parttern occurence start
#'   within a time-series \code{x},
#'   \item at its 2nd element: a vector length of localized parttern occurence
#'   within a time-series \code{x}.
#' }
#'
#' @noRd
#'
finetune_maxima <- function(s.TMP,
                            tau.TMP,
                            nbh.wing,
                            x.Fitted,
                            finetune.maxima.x,
                            template.vl.min,
                            template.vl.max){

  ## Define tau1: pattern occurence start index within a time-series \code{x}
  ## (as preliminarily identified)
  ## Define tau2: pattern occurence end index within a time-series \code{x}
  ## (as preliminarily identified)
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

  ## Identify a pair of points in the two neighbourhods
  ## which corresponds to maxima of `finetune.maxima.x` within egligible indices
  tau1.nbh.x  <- finetune.maxima.x[tau1.nbh]
  tau2.nbh.x  <- finetune.maxima.x[tau2.nbh]
  x.mat       <- outer(tau2.nbh.x, tau1.nbh.x, FUN = "+")
  x.mat.VALID <- x.mat * tau12.mat.VALID
  which.out   <- which(x.mat.VALID == max(x.mat.VALID), arr.ind = TRUE)[1,]

  ## Define "tuned" start and end index point of identified pattern occurence
  ## within a time-series \code{x}
  tau.NEW    <- tau1.nbh[which.out[2]]
  s.NEW      <- tau2.nbh[which.out[1]] - tau.NEW + 1

  return(c(tau.NEW, s.NEW))
}






#' Perform maximization-tunning procedure for ADEPT
#'
#' Perform two-step maximization-tunning procedure to segment pattern occurrences
#' from time-series \code{x}.
#'
#' @details
#' The first step of the procedure consists of maximization of the
#' covariance between the scaled empirical pattern(s) and time-series \code{x}.
#' This provides a good idea about the where the pattern occurrence is localized,
#' but it can miss the exact location by fractions of a second.
#'
#' The second step (optional) is designed to tune the stride segmentation to more match
#' the beginning and ends of the pattern. Currently, only maxima-detection
#' tuning procedure is implemented. Maxima-detection tuning procedure
#' identifies  local maximum of (possibly smoothed) time-series \code{x}
#' in the neighbourhood of preliminary locations of parttern
#' occurence start and end found in step 1.
#'
#' The arguments related to the second (optional) tunning step are:
#' \code{finetune}, \code{finetune.maxima.x}, \code{finetune.maxima.nbh.vl}
#' (all default to \code{NULL}).
#'
#' @param x A numeric vector. Time-series from which we intend to segment pattern occurrences.
#' @param template.vl A numeric vector. A grid of vector lengths that were used for
#' scaling the pattern.
#' @param similarity.mat A numeric matrix with similarity values. A number of matrix columns
#' corresponds to a vector length of time-series \code{x}. A number of matrix rows
#' corresponds to a number of different pattern scale values considered (equivalently:
#' length of \code{template.scaled} list). Each matrix row consists
#' of a vector of similarity statistic (correlation, covariance etc.)
#' between \code{x} (or, possibly, its smoothed version) and a pattern rescaled to matrix row-specific
#' scale parameter; precisely, it is a vector with the highest similarity
#' value corresponding to a particular time point of time-series \code{x}, computed out of
#' possibly multiple patterns.
#' @param similarity.measure.thresh A numeric scalar. Defines threshold of minimal similarity
#' value between time-series \code{x} (or, possibly, its smoothed version) and scaled versions of pattern
#' below which we no longer identify a pattern occurrence.
#' @param finetune A string. Defines type of finetuning procedure empolyed in
#' the procedure. Defaults to \code{NULL}. Currently supported values:
#' \itemize{
#'   \item "maxima" - denotes procedure which tunes preliminary locations of parttern occurence start and end so
#'   as they correspond to local maxima of time-series \code{x} (or, possibly, its smoothed version),
#'   \code{finetune.maxima.x}) found within two respective neighbourhoods of length,
#'   \code{finetune.maxima.nbh.vl}, centered at preliminary locations of parttern occurence start and end.
#' }
#'
#' @param finetune.maxima.x A numeric vector; time-series \code{x} or, possibly, its smoothed version
#' used in "maxima" finetuning procedure.
#' @param finetune.maxima.nbh.vl An integer scalar; denotes length of two respective neighbourhoods,
#'  centered at preliminary locations of parttern occurence start and end, within which the tuned locations of
#'  parttern occurence start and end are searched for.
#'
#' @return A \code{data.frame} object with the segmentation results. Each row
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
#' @noRd
#'
maxAndTune <- function(x,
                       template.vl,
                       similarity.mat,
                       similarity.measure.thresh,
                       finetune = NULL,
                       finetune.maxima.x = NULL,
                       finetune.maxima.nbh.vl = NULL){

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
      finetune.out <- finetune_maxima(s.TMP,
                                      tau.TMP,
                                      nbh.wing,
                                      x.Fitted,
                                      finetune.maxima.x,
                                      template.vl.min,
                                      template.vl.max)
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

