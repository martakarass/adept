


#' ADEPT Similarity Matrix Computation
#'
#' Compute ADEPT similarity matrix between time-series \code{x} windows and a collection
#' of scaled versions of pattern templates.
#'
#' @param x A numerical vector. A time-series \code{x}.
#' @param template.scaled A list of lists of numeric vectors.  Each element of
#' \code{template.scaled} list is a list of pattern templates scaled according
#' to a particular scale parameter (that is, scaled to a particular vector length).
#' Number of elements in the \code{template.scaled} list corresponds to the number
#'  of unique scale values considered in the method. See: \code{scaleTemplate {adept}}
#'  for computation of \code{template.scaled} list.
#' @param similarity.measure A character scalar. Defines a statistic
#' used in similarity matrix computation. One of the following:
#' \itemize{
#'   \item "cov" - for covariance,
#'   \item "cor" - for correlation.
#' }
#'
#' @return A numeric matrix. The matrix entries are similarity values between
#' time-series \code{x} windows and a collection
#' of scaled versions of pattern templates. Each matrix row consists
#' of a vector of similarity statistic between \code{x} windows and a pattern
#' rescaled to a particular scale parameter (that is, scaled to a particular vector length);
#' precisely, for a scale-specific matrix row vector, each vector element is a maximum
#' similarity value selected out of all distinct pattern templates considered.
#'
#' @seealso \code{scaleTemplate {adept}}
#'
#' @export
#' @import runstats
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
#' out <- adeptSimilarity(x, template.scaled, "cov")
#'
#' \dontrun{
#' ## Visualize
#' par(mfrow = c(1,1))
#' image(t(out),
#'       main = "ADEPT similarity matrix\nfor time-series x and scaled versions of a pattern",
#'       xlab = "Time-series x index",
#'       ylab = "Pattern scale (pattern vector length)",
#'       xaxt = "n", yaxt = "n")
#' xaxis <- c(1, seq(1000, length(x0), by = 1000))
#' yaxis <- template.vl
#' axis(1, at = xaxis/max(xaxis), labels = xaxis)
#' axis(2, at = (yaxis - min(yaxis))/(max(yaxis) - min(yaxis)), labels = yaxis)
#' }
#'
adeptSimilarity <- function(x,
                            template.scaled,
                            similarity.measure){

  runstat.func <- switch(similarity.measure,
                         "cov" = RunningCov,
                         "cor" = RunningCor)

  ## Outer lapply: iterate over pattern scales considered;
  ## each lapply iteration fills one row of the output similarity matrix.
  similarity.list <- lapply(template.scaled, function(template.scaled.i){

    ## Inner lapply: iterate over, possibly, multiple patterns;
    ## each lapply iteration returns a vector whose each element corresponds
    ## to the highest value of similarity between signal \code{x} and
    ## a short pattern
    ## at a time point corresponding to this vector's element.
    runstat.func.out0 <- lapply(template.scaled.i, function(template.scaled.ik){
      do.call(runstat.func, list(x = x, y = template.scaled.ik))
    })
    do.call(pmax, runstat.func.out0)

  })

  ## rbind list elements (which are vectors) into a matrix
  similarity.mat <- do.call(rbind, similarity.list)
  return(similarity.mat)

}








