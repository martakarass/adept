


#' ADEPT Similarity Matrix Computation
#'
#' Compute ADEPT similarity matrix between a time-series \code{x} and a collection
#' of scaled templates. Also compute matrix of pattern templates yielding the
#' highest similarity between a time-series \code{x} and a collection
#' of scaled pattern templates.
#'
#' @param x A numeric vector. A time-series \code{x}.
#' @param template.scaled A list of lists of numeric vectors, as returned by
#' \code{scaleTemplate}.  Each element of
#' \code{template.scaled}
#' is a list of templates interpolated to a particular vector length.
#'  Number of elements in the \code{template.scaled}  corresponds to the
#'  number of unique template length values used in segmentation.
#' @param similarity.measure A character scalar. Statistic
#' used in similarity matrix computation; one of the following:
#' \itemize{
#'   \item "cov" - for covariance,
#'   \item "cor" - for correlation.
#' }
#'
#' @return A list of two numeric matrices.
#'
#' "similarity" contains values of similarity between a time-series \code{x}
#' and scaled templates.
#' \itemize{
#'   \item Number of rows equals \code{template.scaled} length,
#'   number of columns equals \code{x}  length.
#'   \item A particular matrix row consists of similarity statistic
#'   between \code{x} and a template rescaled to a particular vector length.
#'   Precisely, each row's element is a maximum out of similarity values
#'  computed for each distinct template used in segmentation.
#' }
#'
#' "idx" represents number of pattern template yielding the highest similarity
#' between a time-series \code{x} and a collection of scaled pattern templates.
#' Precisely, the number is the order in which particular pattern template was
#' provided in the \code{template} list in \code{segmentPattern}.
#'
#' @seealso \code{scaleTemplate {adept}}
#'
#' @export
#'
#' @examples
#' ## Simulate data
#' par(mfrow = c(1,1))
#' x0 <- sin(seq(0, 2 * pi * 100, length.out = 10000))
#' x  <- x0 + rnorm(1000, sd = 0.1)
#' template <- list(x0[1:500])
#' template.vl <- seq(300, 700, by = 50)
#'
#' ## Rescale pattern
#' template.scaled <- scaleTemplate(template, template.vl)
#'
#' ## Compute ADEPT similarity matrix
#' out <- similarityMatrix(x, template.scaled, "cov")$similarity
#'
#' ## Visualize
#' par(mfrow = c(1,1))
#' image(t(out),
#'       main = "ADEPT similarity matrix\nfor time-series x and scaled versions of pattern templates",
#'       xlab = "Time index",
#'       ylab = "Pattern vector length",
#'       xaxt = "n", yaxt = "n")
#' xaxis <- c(1, seq(1000, length(x0), by = 1000))
#' yaxis <- template.vl
#' axis(1, at = xaxis/max(xaxis), labels = xaxis)
#' axis(2, at = (yaxis - min(yaxis))/(max(yaxis) - min(yaxis)), labels = yaxis)
#'
similarityMatrix <- function(x,
                             template.scaled,
                             similarity.measure) {
  if (similarity.measure == "cov") {
    similarity.list <- get_cov_similarities(x, template.scaled)
  } else if (similarity.measure == "cor") {
    similarity.list <- get_cor_similarities(x, template.scaled)
  } else {
    stop("Only 'cov' and 'cor' measures supported")
  }
  ## rbind list elements (which are vectors) into a matrix
  similarity <- do.call(rbind, lapply(similarity.list, function(i) i[[1]]))
  similarity_idx <- do.call(rbind, lapply(similarity.list, function(i) i[[2]]))
  return(list("similarity" = similarity, "idx" = similarity_idx))
}

#' @noRd
#'
get_cov_similarities <- function(x, template.scaled) {
  lapply(template.scaled, function(template.scaled.i) {
    padding <- rep(NA, length(template.scaled.i[[1]]) - 1)
    sliding.func.out0 <-
      lapply(template.scaled.i, function(template.scaled.ik) {
        slidingCovFast(long = x, short = template.scaled.ik)
      })

    maxes <- pmaxIdxCpp(sliding.func.out0)
    list(c(maxes$pmax, padding),
         c(maxes$idx, padding))
  })
}

#' @noRd
#'
get_cor_similarities <- function(x, template.scaled) {
  lapply(template.scaled, function(template.scaled.i) {
    padding <- rep(NA, length(template.scaled.i[[1]]) - 1)

    first_template <-
      slidingCorStoreSd(long = x, short = template.scaled.i[[1]])
    if (length(template.scaled.i) > 1) {
      sliding.func.out0 <-
        lapply(template.scaled.i[2:length(template.scaled.i)], function(template.scaled.ik) {
          slidingCor(long = x,
                     short = template.scaled.ik,
                     sds = first_template$sds)
        })
      maxes <-
        pmaxIdxCpp(c(list(first_template$cor), sliding.func.out0))
      return(list(c(maxes$pmax, padding),
                  c(maxes$idx, padding)))
    } else {
      return(list(c(first_template$cor, padding),
                  c(rep(
                    1, length(first_template$cor)
                  ), padding)))
    }
  })
}
