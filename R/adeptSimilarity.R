


#' Title
#'
#' @param x
#' @param template.scaled
#' @param similarity.measure
#'
#' @return
#'
#' @import runstats
#'
#' @export
#'
#' @examples
#'
adeptSimilarity <- function(x,
                            template.scaled,
                            similarity.measure){

  runstat.func <- switch(similarity.measure,
                         "cov" = RunningCov,
                         "cor" = RunningCor)

  ## Outer lapply: iterate over template scales considered:
  ## each lapply iteration "fills one row" of final similarity matrix
  similarity.list <- lapply(template.scaled, function(template.scaled.i){

    ## Outer lapply: iterate over distinct templates considered (fixed template scale)
    runstat.func.out0 <- lapply(template.scaled.i, function(template.scaled.ik){
      do.call(runstat.func, list(x = x, y = template.scaled.ik))
    })

    ## Compute and store highest result over possibly multiple templates
    do.call(pmax, runstat.func.out0)
  })

  ## rbind into a matrix
  similarity.mat <- do.call(rbind, similarity.list)
  return(similarity.mat)

}
