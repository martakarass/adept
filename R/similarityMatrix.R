


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
                            similarity.measure){

  sliding.func <- switch(similarity.measure,
                         "cov" = slidingCovFast,
                         "cor" = slidingCorFast)

  ## Outer lapply: iterate over pattern scales considered;
  ## each lapply iteration fills one row of the output similarity matrix.
  similarity.list <- lapply(template.scaled, function(template.scaled.i){

    ## Inner lapply: iterate over, possibly, multiple patterns;
    ## each lapply iteration returns a vector whose each element corresponds
    ## to the highest value of similarity between signal \code{x} and
    ## a short pattern
    ## at a time point corresponding to this vector's element.
    sliding.func.out0 <- lapply(template.scaled.i, function(template.scaled.ik){
      do.call(sliding.func, list(long = x, short = template.scaled.ik))
    })
    maxes <- pmaxIdxCpp(sliding.func.out0)
    padding <- rep(NA, length(template.scaled.i[[1]]) - 1)
    list(c(maxes$pmax, padding),
         c(maxes$idx, padding))
  })

  similarity <- do.call(rbind, lapply(similarity.list, function(i) i[[1]]))
  similarity.idx <- do.call(rbind, lapply(similarity.list, function(i) i[[2]]))
  return(list("similarity" = similarity, "idx" = similarity.idx))
}
