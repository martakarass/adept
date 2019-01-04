


#' Compute ADEPT similarity matrix
#'
#' Compute ADEPT similarity matrix between time-series \code{x} and a collection
#' of rescaled versions of empirical pattern(s).
#'
#' @param x A numeric vector; time-series \code{x}.
#' @param template.scaled A list of rescaled versions of empirical pattern(s).
#' Number of elements in the list corresponds to number of unique pattern scale
#' parameters considered. Each element of \code{template.scaled} is a list
#' itself and contains numeric vector(s) where each vector corresponds to
#' one of possibly multiple distinct patterns considered.
#' @param similarity.measure Character value; one of the following:
#' \itemize{
#'   \item "cov",
#'   \item "cor".
#' }
#'
#' @return A numeric matrix with similarity values; number of columns
#' corresponds to vector length of time-series \code{x}. Number of rows
#' corresponds different pattern scale values considered (equivalently:
#' length of \code{template.scaled} list). Each row of the matrix consists
#' of vector of running similarity statistic (correlation, covariance etc.)
#' between  \code{x} and a pattern template(s) rescaled to matrix row-specific
#' particular scale parameter. Precisely, a vector with the highest similarity
#' result at each time point, over possibly multiple pattern templates, is
#' returned.
#'
#' @import runstats
#'
#' @noRd
#'
adeptSimilarity <- function(x,
                            template.scaled,
                            similarity.measure){

  runstat.func <- switch(similarity.measure,
                         "cov" = RunningCov,
                         "cor" = RunningCor)

  ## Outer lapply: iterate over template scales considered:
  ## each lapply iteration fills one row of an output similarity matrix
  similarity.list <- lapply(template.scaled, function(template.scaled.i){

    ## Inner lapply: iterate over distinct pattern templates considered
    ## (fixed template scale); a vector with the highest similarity result
    ## at each time point, over possibly multiple pattern templates, is returned
    runstat.func.out0 <- lapply(template.scaled.i, function(template.scaled.ik){
      do.call(runstat.func, list(x = x, y = template.scaled.ik))
    })
    do.call(pmax, runstat.func.out0)

  })

  ## rbind list elements (which are vectors) into a matrix
  similarity.mat <- do.call(rbind, similarity.list)
  return(similarity.mat)

}
